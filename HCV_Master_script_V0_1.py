# HCV_Master_script_V0_1.py
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path
import sys
import os
import glob
import re
import datetime
from Bio import SeqIO
from collections import defaultdict, deque
import gzip
import numpy as np

def run_script(script_path: Path, args_list: list, capture_stdout=True, cwd=None) -> tuple[int, str | None]:
    """Runs a Python script using subprocess and returns status code and stdout."""
    command = [sys.executable, str(script_path)] + args_list
    print(f"Running command: {' '.join(command)}", file=sys.stderr)
    try:
        process = subprocess.run(command, check=False, text=True, capture_output=True, cwd=cwd)
        if process.returncode != 0:
            print(f"Error: Script {script_path.name} failed with exit code {process.returncode}", file=sys.stderr)
            if process.stderr:
                print(f"Stderr:\n{process.stderr}", file=sys.stderr)
            if process.stdout and not capture_stdout:
                 print(f"Stdout:\n{process.stdout}", file=sys.stderr)
        output_to_return = process.stdout if capture_stdout else (process.stderr if process.returncode != 0 else None)
        return process.returncode, output_to_return
    except FileNotFoundError:
        print(f"Error: Script not found: {script_path}", file=sys.stderr)
        return 1, None
    except Exception as e:
        print(f"Error running script {script_path}: {e}", file=sys.stderr)
        return 1, None

def find_clusters(linked_pairs: set) -> list[set]:
    """Finds connected components (clusters) in a graph defined by linked pairs."""
    adj = defaultdict(set)
    all_samples = set()
    for u, v in linked_pairs:
        adj[u].add(v)
        adj[v].add(u)
        all_samples.add(u)
        all_samples.add(v)

    visited = set()
    clusters = []
    for sample in sorted(list(all_samples)): # Process samples consistently
        if sample not in visited:
            current_cluster = set()
            q = deque([sample])
            visited.add(sample)
            while q:
                u = q.popleft()
                current_cluster.add(u)
                for v in adj[u]:
                    if v not in visited:
                        visited.add(v)
                        q.append(v)
            if current_cluster:
                clusters.append(current_cluster)
    return clusters

def load_existing_clusters(cluster_file_path: Path) -> tuple[dict[str, set[str]], dict[str, str], int, set[tuple[str, str]]]:
    """
    Loads existing clusters from a file.
    Returns:
        - existing_clusters_map: Dict mapping cluster_name -> set(samples)
        - sample_to_cluster_map: Dict mapping sample -> cluster_name
        - max_cluster_num: The highest cluster number found (e.g., 2 for Cluster_2)
        - existing_links: Set of tuples representing links implied by cluster membership
    """
    existing_clusters_map = {}
    sample_to_cluster_map = {}
    existing_links = set()
    max_cluster_num = 0
    if cluster_file_path.is_file():
        print(f"Loading existing clusters from: {cluster_file_path}", file=sys.stderr)
        try:
            with open(cluster_file_path, "r") as f:
                for line in f:
                    if line.strip() and not line.startswith("#"):
                        parts = line.strip().split('\t')
                        if len(parts) == 2:
                            cluster_name, samples_str = parts
                            try:
                                cluster_num = int(cluster_name.split('_')[-1])
                                max_cluster_num = max(max_cluster_num, cluster_num)
                            except (IndexError, ValueError):
                                print(f"Warning: Could not parse cluster number from '{cluster_name}' in {cluster_file_path}", file=sys.stderr)
                                continue # Skip this line

                            samples = [s.strip() for s in samples_str.split(',') if s.strip()]
                            cluster_sample_set = set(samples)
                            existing_clusters_map[cluster_name] = cluster_sample_set
                            for sample in samples:
                                sample_to_cluster_map[sample] = cluster_name
                            # Add links between all samples in the existing cluster
                            if len(samples) >= 2:
                                for i in range(len(samples)):
                                    for j in range(i + 1, len(samples)):
                                        existing_links.add(tuple(sorted((samples[i], samples[j]))))
                        else:
                             print(f"Warning: Skipping malformed line in {cluster_file_path}: {line.strip()}", file=sys.stderr)
        except IOError as e:
            print(f"Warning: Could not read existing cluster file {cluster_file_path}: {e}", file=sys.stderr)
    # next_cluster_id logic moved to main clustering section
    print(f"Loaded {len(existing_clusters_map)} existing clusters. Max cluster number: {max_cluster_num}. Implied links: {len(existing_links)}.", file=sys.stderr)
    return existing_clusters_map, sample_to_cluster_map, max_cluster_num, existing_links


def main():
    parser = argparse.ArgumentParser(description="Master script for HCV transmission pipeline.")
    parser.add_argument("reads_dir", type=Path, help="Directory containing input FASTQ files (e.g., .../reads/).")
    parser.add_argument("base_dir", type=Path, help="Base directory for the pipeline scripts and data subdirectories (Reports, HCV_Kmers, etc.).")
    parser.add_argument("--keep_tmp", action='store_true', help="Keep temporary run directories.")
    parser.add_argument("--keep_unmerged", action='store_true', help="Keep unmerged reads from the combine step.") # New argument
    parser.add_argument("--combine_script", type=str, default="HCV_combine_reads_V0_1.py", help="Name of the combine reads script.")
    parser.add_argument("--distancer_script", type=str, default="HCV_kmers_distancer_V0_1.py", help="Name of the kmer distancer script.")
    parser.add_argument("--transmission_script", type=str, default="HCV_transmission_test_V0_2.py", help="Name of the transmission test script.")
    parser.add_argument("--r1_pattern", type=str, default="*_R1_001.fastq.gz", help="Glob pattern for R1 files.")
    parser.add_argument("--kmer_overlap_threshold", type=float, default=0.05, help="Threshold for k-mer overlap ratio to trigger transmission test.")
    parser.add_argument("--snp_dists_path", type=str, default="snp-dists", help="Path to the snp-dists executable.")
    parser.add_argument("--rscript_path", type=str, default="Rscript", help="Path to the Rscript executable.")
    args = parser.parse_args()

    # --- Validate paths and setup directories ---
    if not args.reads_dir.is_dir():
        print(f"Error: Reads directory not found: {args.reads_dir}", file=sys.stderr)
        sys.exit(1)
    if not args.base_dir.is_dir():
        print(f"Error: Base directory not found: {args.base_dir}", file=sys.stderr)
        sys.exit(1)

    reports_dir = args.base_dir / "Reports"
    kmer_dir = args.base_dir / "HCV_Kmers"
    fasta_dir = args.base_dir / "HCV_fasta"
    combined_reads_dir = args.base_dir / "combined_reads" # combine script creates this

    reports_dir.mkdir(parents=True, exist_ok=True)

    combine_script_path = args.base_dir / args.combine_script
    distancer_script_path = args.base_dir / args.distancer_script
    transmission_script_path = args.base_dir / args.transmission_script # Path to the next script

    if not combine_script_path.is_file():
         print(f"Error: Combine script not found: {combine_script_path}", file=sys.stderr)
         sys.exit(1)
    if not distancer_script_path.is_file():
         print(f"Error: Distancer script not found: {distancer_script_path}", file=sys.stderr)
         sys.exit(1)
    # Transmission script existence check moved to where it's needed

    processed_samples = []
    files_to_send = [] # For potential emailing
    master_log_file = reports_dir / "master_pipeline.log" # Central log file
    run_report_lines = [] # Store lines for the final run report

    # --- Load Existing Clusters ---
    cluster_file_path = reports_dir / "clusters.txt"
    # Load existing cluster state
    existing_clusters_map, initial_sample_to_cluster, max_cluster_num, existing_links = load_existing_clusters(cluster_file_path)
    # Combine existing links with newly found links for downstream processing if needed, but clustering logic below uses only *new* links to modify state
    all_linked_pairs = existing_links.copy() # Start with links implied by existing clusters
    newly_linked_pairs_this_run = set() # Track only links found *in this specific run*

    # --- Step 1: Combine Reads ---
    print("--- Running Step 1: Combine Reads ---", file=sys.stderr)
    r1_files = list(args.reads_dir.glob(args.r1_pattern))

    if not r1_files:
        print(f"Warning: No files found matching pattern '{args.r1_pattern}' in {args.reads_dir}", file=sys.stderr)

    # Dictionary to store the main temp_run_dir for each processed sample
    sample_temp_dirs = {}

    for r1_file in r1_files:
        if "Undeter" in r1_file.name or "undeter" in r1_file.name:
            print(f"Skipping undetermined file: {r1_file.name}", file=sys.stderr)
            continue

        # Determine R2
        r1_filename = r1_file.name
        r2_file = None
        try:
            if "_R1_" in r1_filename: r2_filename = r1_filename.replace("_R1_", "_R2_", 1)
            elif "_R1." in r1_filename: r2_filename = r1_filename.replace("_R1.", "_R2.", 1)
            elif "_1." in r1_filename: r2_filename = r1_filename.replace("_1.", "_2.", 1)
            elif "_1_" in r1_filename: r2_filename = r1_filename.replace("_1_", "_2_", 1)
            else: raise ValueError("Cannot determine R2 pattern")
            r2_file_check = r1_file.parent / r2_filename
            if not r2_file_check.is_file():
                 print(f"Warning: R2 file not found for {r1_file.name}, skipping.", file=sys.stderr)
                 continue
            r2_file = r2_file_check
            prefix = r1_filename.split('_')[0]
            if not prefix: prefix = r1_file.stem
        except Exception as e:
            print(f"Warning: Could not determine R2 or prefix for {r1_file.name}: {e}. Skipping.", file=sys.stderr)
            continue

        print(f"\nProcessing sample: {prefix} ({r1_file.name})", file=sys.stderr)
        report_file = reports_dir / f"{prefix}_report.out"
        files_to_send.append(report_file)
        run_report_lines.append(f"--- Sample: {prefix} ---") # Add sample to run report

        temp_run_dir = None
        try:
            temp_run_dir_path = tempfile.mkdtemp(prefix=f'run_{prefix}_', dir=args.base_dir)
            temp_run_dir = Path(temp_run_dir_path)
            sample_temp_dirs[prefix] = temp_run_dir
            print(f"\nCreated temporary directory: {temp_run_dir}\n", file=sys.stderr)

            with open(report_file, 'w') as f_out:
                f_out.write(f"Processing sample {prefix}:\n")
                f_out.write(f"R1 file: {r1_file}\n")
                f_out.write(f"R2 file: {r2_file}\n")
                f_out.write(f"Timestamp: {datetime.datetime.now()}\n\n")
                f_out.write("--- Combine Reads Output ---\n")

            combine_args = [str(args.base_dir), str(r1_file), str(temp_run_dir)]
            # Pass down the keep_unmerged flag if set
            if args.keep_unmerged:
                combine_args.append("--keep_unmerged")

            return_code, stdout = run_script(combine_script_path, combine_args, capture_stdout=True)

            with open(report_file, 'a') as f_out:
                if return_code == 0 and stdout:
                    f_out.write(stdout)
                    processed_samples.append(prefix)
                else:
                    f_out.write(f"\nERROR: Combine script failed with exit code {return_code}.\n")
                    print(f"\nERROR: Combine script failed for {prefix}. See report file.\n", file=sys.stderr)
                    run_report_lines.append(f"  Status: Combine Reads FAILED (Code: {return_code})")
                f_out.write("\n--- End Combine Reads Output ---\n\n")

        except Exception as e:
            print(f"Error during combine step for {prefix}: {e}", file=sys.stderr)
            run_report_lines.append(f"  Status: Combine Reads FAILED (Exception: {e})")
            if report_file.exists():
                 with open(report_file, 'a') as f_out: f_out.write(f"\nFATAL ERROR during combine step: {e}\n")
        # finally block moved to end of Step 2

    
    # --- Step 2: Kmer Distancing and SNP Linking ---
    
    print("\n--- Running Step 2: Kmer Distancing and SNP Linking ---", file=sys.stderr)
    if not processed_samples:
        print("No samples successfully processed in Step 1. Skipping Step 2.", file=sys.stderr)
    else:
        print(f"Analyzing {len(processed_samples)} samples: {', '.join(processed_samples)}", file=sys.stderr)

    sample_link_info = defaultdict(list) # Store linked sample and min distance for the run report: {sample: [(linked_sample, min_dist), ...]}

    for sample_prefix in processed_samples:
        print(f"\nAnalyzing k-mers for sample: {sample_prefix}", file=sys.stderr)
        report_file = reports_dir / f"{sample_prefix}_report.out"
        kmer_file = kmer_dir / f"{sample_prefix}.kmers.gz"
        fasta_file = fasta_dir / f"{sample_prefix}.fasta.gz"
        temp_run_dir = sample_temp_dirs.get(sample_prefix)

        if not temp_run_dir:
             print(f"Error: Could not find temporary directory for sample {sample_prefix}. Skipping.", file=sys.stderr)
             run_report_lines.append(f"  Status: SKIPPED (Missing temp dir)")
             continue

        try:
            with open(report_file, 'a') as f_out:
                f_out.write("\n--- Kmer Distance Analysis ---\n")
                if not kmer_file.is_file():
                    f_out.write(f"ERROR: Kmer file not found: {kmer_file}\nSkipping distance analysis.\n")
                    print(f"Warning: Kmer file not found for {sample_prefix}, skipping distance analysis.", file=sys.stderr)
                    run_report_lines.append(f"  Status: Kmer Analysis SKIPPED (No kmer file)")
                    continue

                distancer_args = [str(args.base_dir), sample_prefix, f"--overlap_threshold={args.kmer_overlap_threshold}"]
                return_code, stdout = run_script(distancer_script_path, distancer_args, capture_stdout=True)

                if return_code != 0 or not stdout:
                    f_out.write(f"ERROR: Kmer distancer script failed with exit code {return_code}.\n")
                    print(f"ERROR: Kmer distancer failed for {sample_prefix}.", file=sys.stderr)
                    run_report_lines.append(f"  Status: Kmer Analysis FAILED (Code: {return_code})")
                    continue

                f_out.write(stdout)
                f_out.write("\n--- End Kmer Distance Analysis ---\n")

                transmission_candidates = []
                potential_links = 0
                lines = stdout.strip().split('\n')
                if len(lines) > 1:
                    for line in lines[1:-1]:
                        parts = line.strip().split('\t')
                        if len(parts) == 6:
                            try:
                                ratio = float(parts[5].rstrip('*'))
                                if ratio >= args.kmer_overlap_threshold:
                                    potential_links += 1
                                    transmission_candidates.append(parts[1])
                            except (ValueError, IndexError):
                                print(f"Warning: Could not parse kmer distance line: {line}", file=sys.stderr)
                                f_out.write(f"\nWarning: Could not parse kmer distance line: {line}\n")

                f_out.write(f"\n--- SNP Distance Check (Threshold <= 9) ---\n")
                if not transmission_candidates:
                    f_out.write(f"No potential transmission candidates found based on k-mer overlap >= {args.kmer_overlap_threshold}.\n")
                    run_report_lines.append(f"  Links found: No (No k-mer candidates)")
                else:
                    f_out.write(f"{potential_links} k-mer candidate link(s) detected for sample {sample_prefix}.\n")
                    found_any_snp_link_for_sample = False # Track if any link was found for this sample_prefix

                    if not Path(args.snp_dists_path).is_file() and shutil.which(args.snp_dists_path) is None:
                         f_out.write(f"ERROR: snp-dists executable not found at '{args.snp_dists_path}'. Cannot perform SNP checks.\n")
                         print(f"ERROR: snp-dists executable not found: {args.snp_dists_path}", file=sys.stderr)
                         run_report_lines.append(f"  Links found: Error (snp-dists not found)")
                         continue

                    if not fasta_file.is_file():
                         f_out.write(f"ERROR: Query FASTA file not found: {fasta_file}. Cannot perform SNP checks.\n")
                         print(f"Warning: Query FASTA file not found for {sample_prefix}, skipping SNP checks.", file=sys.stderr)
                         run_report_lines.append(f"  Links found: Error (Query FASTA not found)")
                         continue

                    for subject_prefix in transmission_candidates:
                        subject_fasta_file = fasta_dir / f"{subject_prefix}.fasta.gz"
                        if not subject_fasta_file.is_file():
                            f_out.write(f"ERROR: Subject FASTA file not found: {subject_fasta_file}. Skipping SNP check for this pair.\n")
                            print(f"Warning: Subject FASTA file not found for {subject_prefix}, skipping pair {sample_prefix} vs {subject_prefix}.", file=sys.stderr)
                            continue

                        # Define paths for intermediate files within the main temp_run_dir
                        pair_combined_fasta = temp_run_dir / f"pair_{sample_prefix}_{subject_prefix}_combined.fasta"
                        pair_aligned_fasta = temp_run_dir / f"pair_{sample_prefix}_{subject_prefix}_aligned.fasta"
                        pair_unique_headers_fasta = temp_run_dir / f"pair_{sample_prefix}_{subject_prefix}_unique_headers.fasta"
                        pair_snp_list_file = temp_run_dir / f"snp_distances_{sample_prefix}_{subject_prefix}.tsv" # Output file for snp-dists -m

                        try:
                            pair_headers_s1 = set()
                            pair_headers_s2 = set()
                            with gzip.open(fasta_file, "rt") as handle:
                                for record in SeqIO.parse(handle, "fasta"): pair_headers_s1.add(record.id.replace('_R_', ''))
                            with gzip.open(subject_fasta_file, "rt") as handle:
                                for record in SeqIO.parse(handle, "fasta"): pair_headers_s2.add(record.id.replace('_R_', ''))

                            with open(pair_combined_fasta, "wb") as f_cat:
                                for f_gz in [fasta_file, subject_fasta_file]:
                                    with gzip.open(f_gz, "rb") as f_in: shutil.copyfileobj(f_in, f_cat)

                            mafft_docker_cmd = ["docker", "run", "--rm", f"-v", f"{temp_run_dir.resolve()}:/data", "pegi3s/mafft", "mafft", "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder", f"/data/{pair_combined_fasta.name}"]
                            print(f"Running MAFFT for pair {sample_prefix} vs {subject_prefix}...", file=sys.stderr)
                            mafft_ret = 1
                            try:
                                with open(pair_aligned_fasta, 'w') as f_out_aln:
                                    proc_mafft = subprocess.run(mafft_docker_cmd, cwd=temp_run_dir, text=True, capture_output=True, check=True)
                                    f_out_aln.write(proc_mafft.stdout)
                                mafft_ret = 0
                            except Exception as e:
                                print(f"MAFFT failed for pair {sample_prefix} vs {subject_prefix}: {e}", file=sys.stderr)
                                f_out.write(f"ERROR: MAFFT failed for pair {sample_prefix} vs {subject_prefix}. Cannot check SNPs.\n")
                                continue

                            if mafft_ret != 0 or not pair_aligned_fasta.is_file() or pair_aligned_fasta.stat().st_size == 0:
                                print(f"MAFFT failed or produced empty alignment for pair {sample_prefix} vs {subject_prefix}. Skipping SNP check.", file=sys.stderr)
                                f_out.write(f"ERROR: MAFFT alignment failed for pair {sample_prefix} vs {subject_prefix}. Cannot check SNPs.\n")
                                continue

                            pair_unique_id_to_sample = {}
                            with open(pair_aligned_fasta, "r") as f_in_aln, open(pair_unique_headers_fasta, "w") as f_out_unique:
                                for record in SeqIO.parse(f_in_aln, "fasta"):
                                    original_id_cleaned = record.id.replace('_R_', '')
                                    current_sample = None
                                    if original_id_cleaned in pair_headers_s1: current_sample = sample_prefix
                                    elif original_id_cleaned in pair_headers_s2: current_sample = subject_prefix
                                    if current_sample: record.id = f"{current_sample}_{record.id}"; pair_unique_id_to_sample[record.id] = current_sample
                                    else: record.id = f"UNKNOWN_{record.id}"; pair_unique_id_to_sample[record.id] = "UNKNOWN"
                                    record.description = ""; SeqIO.write(record, f_out_unique, "fasta")

                            snp_dists_cmd_str = f"{args.snp_dists_path} -m {pair_unique_headers_fasta.name} > {pair_snp_list_file.name}"
                            print(f"Running snp-dists for pair {sample_prefix} vs {subject_prefix}...", file=sys.stderr)
                            snp_dist_ret = 1
                            try:
                                proc_snp = subprocess.run(snp_dists_cmd_str, shell=True, cwd=temp_run_dir, text=True, capture_output=True, check=False)
                                snp_dist_ret = proc_snp.returncode
                                if snp_dist_ret != 0:
                                    print(f"snp-dists failed for pair {sample_prefix} vs {subject_prefix} (Exit Code {snp_dist_ret}).", file=sys.stderr)
                                    if proc_snp.stderr: print(f"Stderr:\n{proc_snp.stderr}", file=sys.stderr)
                                    f_out.write(f"ERROR: snp-dists failed for pair {sample_prefix} vs {subject_prefix}.\n")
                            except Exception as e:
                                print(f"Error running snp-dists for pair {sample_prefix} vs {subject_prefix}: {e}", file=sys.stderr)
                                f_out.write(f"ERROR: Exception running snp-dists for pair {sample_prefix} vs {subject_prefix}: {e}\n")

                            snp_link_found = False
                            min_snp_dist_for_pair = float('inf')
                            inter_pair_count = 0
                            le_10_count = 0
                            if snp_dist_ret == 0 and pair_snp_list_file.is_file() and pair_snp_list_file.stat().st_size > 0:
                                print(f"Checking SNP distances in {pair_snp_list_file.name} (3-column format) for link <= 9...", file=sys.stderr)
                                try:
                                    with open(pair_snp_list_file, "r") as f_pairs:
                                        header = f_pairs.readline()
                                        if not header.startswith("Sample1\tSample2"): f_pairs.seek(0)
                                        for line_num, line in enumerate(f_pairs, start=1):
                                            parts = line.strip().split('\t')
                                            if len(parts) != 3: continue
                                            seq1_id, seq2_id, dist_str = parts
                                            sample1 = pair_unique_id_to_sample.get(seq1_id)
                                            sample2 = pair_unique_id_to_sample.get(seq2_id)
                                            is_target_pair = (sample1 == sample_prefix and sample2 == subject_prefix) or (sample1 == subject_prefix and sample2 == sample_prefix)
                                            if is_target_pair:
                                                inter_pair_count += 1
                                                try:
                                                    dist = int(dist_str)
                                                    min_snp_dist_for_pair = min(min_snp_dist_for_pair, dist)
                                                    if dist <= 9:
                                                        le_10_count += 1
                                                        # Record the link and update minimum distance regardless of whether it's the *first* link found <=9
                                                        # This ensures we capture the true minimum if multiple pairs are <= 9
                                                        # print(f"Found SNP link candidate (dist={dist} <= 9) between {sample_prefix} and {subject_prefix}", file=sys.stderr)
                                                        link_tuple = tuple(sorted((sample_prefix, subject_prefix)))
                                                        # Add to *both* sets: all_linked_pairs for potential downstream use,
                                                        # and newly_linked_pairs_this_run for the stateful clustering logic
                                                        all_linked_pairs.add(link_tuple)
                                                        newly_linked_pairs_this_run.add(link_tuple) # Track new links
                                                        snp_link_found = True # Mark that *a* link <= 9 was found for this pair
                                                        # We already update min_snp_dist_for_pair in line 413 (corrected line number)
                                                except ValueError:
                                                    print(f"Warning: Could not parse distance '{dist_str}' as integer in {pair_snp_list_file.name} line {line_num}", file=sys.stderr)
                                except Exception as e:
                                    print(f"Error parsing {pair_snp_list_file.name}: {e}", file=sys.stderr)
                                    f_out.write(f"\nERROR: Could not parse SNP distance list: {e}\n")

                            min_dist_report = str(min_snp_dist_for_pair) if min_snp_dist_for_pair != float('inf') else "N/A"
                            percentage = (le_10_count / inter_pair_count * 100) if inter_pair_count > 0 else 0
                            f_out.write(f"Minimum SNP distance found between {sample_prefix} and {subject_prefix}: {min_dist_report}\n")
                            f_out.write(f"Percentage of pairs <= 9 SNPs: {percentage:.2f}% ({le_10_count}/{inter_pair_count})\n")
                            if snp_link_found:
                                f_out.write(f"** Link confirmed by minimum distance ({min_snp_dist_for_pair} <= 9 SNPs) **\n")
                                # Store the linked sample AND the minimum distance found for this specific pair comparison
                                # Check if this subject is already linked, update distance if lower
                                found_existing = False
                                for i, (existing_subj, existing_dist) in enumerate(sample_link_info[sample_prefix]):
                                    if existing_subj == subject_prefix:
                                        if min_snp_dist_for_pair < existing_dist:
                                            sample_link_info[sample_prefix][i] = (subject_prefix, min_snp_dist_for_pair)
                                        found_existing = True
                                        break
                                if not found_existing:
                                     sample_link_info[sample_prefix].append((subject_prefix, min_snp_dist_for_pair))
                                found_any_snp_link_for_sample = True

                        except Exception as pair_proc_e:
                             print(f"Error during MAFFT/snp-dists processing for pair {sample_prefix} vs {subject_prefix}: {pair_proc_e}", file=sys.stderr)
                             f_out.write(f"ERROR during MAFFT/snp-dists processing for pair: {pair_proc_e}\n")

                        f_out.write("\n") # Add spacing between pairs

                    # Add summary to run report lines *after* checking all candidates for this sample_prefix
                    if not transmission_candidates:
                         run_report_lines.append(f"  Links found: No (No k-mer candidates)")
                    elif not found_any_snp_link_for_sample:
                         run_report_lines.append(f"  Links found: No (SNP distance > 9 or errors)")
                    else:
                         link_summaries = []
                         # Use the updated sample_link_info which now stores (linked_sample, min_dist) tuples
                         if sample_prefix in sample_link_info:
                             for linked_samp, min_dist in sorted(sample_link_info[sample_prefix]): # Sort for consistent output
                                  min_dist_str = str(min_dist) if min_dist != float('inf') else "xx" # Use "xx" as per desired format
                                  link_summaries.append(f"Links found to sample {linked_samp}. Minimum SNP distance is {min_dist_str}.") # Match desired format
                             if link_summaries:
                                  run_report_lines.extend(link_summaries) # Add each link as a separate line
                         else: # Should not happen if found_any_snp_link_for_sample is True, but as fallback
                              run_report_lines.append(f"  Links found: Yes (Details unavailable)") # Fallback message


                f_out.write("\n--- End SNP Distance Check ---\n") # Changed footer

        except Exception as e:
            print(f"Error during Step 2 (kmer/SNP check) for {sample_prefix}: {e}", file=sys.stderr)
            run_report_lines.append(f"  Status: FAILED during Kmer/SNP Check (Exception: {e})")
            if report_file.exists():
                 with open(report_file, 'a') as f_out:
                      f_out.write(f"\nFATAL ERROR during kmer/SNP check analysis: {e}\n")

        # Clean up main temp_run_dir for the sample *after* all pairs have been checked
        finally:
            if temp_run_dir and temp_run_dir.exists() and not args.keep_tmp:
                print(f"Removing temporary directory for sample {sample_prefix}: {temp_run_dir}", file=sys.stderr)
                shutil.rmtree(temp_run_dir, ignore_errors=True)
            elif temp_run_dir and temp_run_dir.exists() and args.keep_tmp:
                print(f"Keeping temporary directory for sample {sample_prefix}: {temp_run_dir}", file=sys.stderr)


    # --- Step 3: Stateful Clustering Update ---
    print("\n--- Running Step 3: Stateful Clustering Update ---", file=sys.stderr)
    # Initialize current state from loaded data
    current_clusters = existing_clusters_map.copy() # {cluster_name: {samples}}
    current_sample_to_cluster = initial_sample_to_cluster.copy() # {sample: cluster_name}
    next_cluster_id_counter = max_cluster_num + 1
    clusters_to_delete = set() # Track names of clusters removed by merging

    # Track changes for reporting
    report_events = defaultdict(list) # {sample: ["event description", ...]}
    merged_cluster_log = {} # {removed_cluster_name: merged_into_cluster_name}
    modified_clusters_this_run = set() # Track clusters created/modified in this run

    print(f"Processing {len(newly_linked_pairs_this_run)} new links found in this run.", file=sys.stderr)
    new_cluster_candidates = set() # Samples involved in new links but not in existing clusters
    new_links_for_new_clusters = set() # Links between samples not initially in clusters

    # --- Process newly found links to update cluster state ---
    for s1, s2 in newly_linked_pairs_this_run:
        c1 = current_sample_to_cluster.get(s1)
        c2 = current_sample_to_cluster.get(s2)

        # Case 1: Link merges two existing, different clusters
        if c1 and c2 and c1 != c2 and c1 not in clusters_to_delete and c2 not in clusters_to_delete:
            num1 = int(c1.split('_')[-1])
            num2 = int(c2.split('_')[-1])
            if num1 < num2:
                keep_c, remove_c = c1, c2
            else:
                keep_c, remove_c = c2, c1

            if remove_c not in clusters_to_delete: # Ensure we haven't already decided to remove this one
                print(f"Merging cluster {remove_c} into {keep_c} due to link between {s1} and {s2}", file=sys.stderr)
                samples_to_move = current_clusters.get(remove_c, set())
                current_clusters[keep_c].update(samples_to_move)
                modified_clusters_this_run.add(keep_c) # Mark target cluster as modified
                for sample in samples_to_move:
                    current_sample_to_cluster[sample] = keep_c # Update mapping
                clusters_to_delete.add(remove_c)
                merged_cluster_log[remove_c] = keep_c
                # Report event for samples involved in the *link* causing the merge
                final_members = sorted(list(current_clusters[keep_c]))
                merge_msg = f"caused merge of {remove_c} into {keep_c} -> {', '.join(final_members)}"
                report_events[s1].append(merge_msg)
                report_events[s2].append(merge_msg)
                # We might want a more global merge report message later

        # Case 2: Link adds a new sample to an existing cluster
        elif c1 and not c2 and c1 not in clusters_to_delete: # s1 in cluster, s2 is new
            print(f"Adding sample {s2} to cluster {c1} due to link with {s1}", file=sys.stderr)
            current_clusters[c1].add(s2)
            current_sample_to_cluster[s2] = c1
            modified_clusters_this_run.add(c1) # Mark cluster as modified
            report_events[s2].append(f"added to existing cluster {c1}")
        elif c2 and not c1 and c2 not in clusters_to_delete: # s2 in cluster, s1 is new
            print(f"Adding sample {s1} to cluster {c2} due to link with {s2}", file=sys.stderr)
            current_clusters[c2].add(s1)
            current_sample_to_cluster[s1] = c2
            modified_clusters_this_run.add(c2) # Mark cluster as modified
            report_events[s1].append(f"added to existing cluster {c2}")

        # Case 3: Link is between samples already in the same cluster (or a cluster marked for deletion)
        elif c1 and c2 and c1 == c2:
             # No structural change, but note for report that they are linked within the cluster
             if c1 not in clusters_to_delete: # Only report if the cluster still exists
                 report_events[s1].append(f"linked within {c1}")
                 report_events[s2].append(f"linked within {c1}")
             pass # Already in the same cluster

        # Case 4: Link is between two samples not currently in any cluster
        elif not c1 and not c2:
            new_cluster_candidates.add(s1)
            new_cluster_candidates.add(s2)
            new_links_for_new_clusters.add(tuple(sorted((s1, s2))))

    # --- Find and create completely new clusters from remaining candidates ---
    if new_links_for_new_clusters:
        print(f"Finding new clusters among {len(new_cluster_candidates)} candidates using {len(new_links_for_new_clusters)} links.", file=sys.stderr)
        newly_formed_clusters = find_clusters(new_links_for_new_clusters)
        for new_cluster_set in newly_formed_clusters:
            new_cluster_name = f"Cluster_{next_cluster_id_counter}"
            print(f"Creating new {new_cluster_name} with samples: {', '.join(sorted(list(new_cluster_set)))}", file=sys.stderr)
            current_clusters[new_cluster_name] = new_cluster_set
            modified_clusters_this_run.add(new_cluster_name) # Mark new cluster
            for sample in new_cluster_set:
                current_sample_to_cluster[sample] = new_cluster_name
                report_events[sample].append(f"formed new cluster {new_cluster_name}")
            next_cluster_id_counter += 1

    # --- Write the final cluster state to file ---
    cluster_file_path = reports_dir / "clusters.txt"
    final_cluster_count = 0
    cluster_info_for_plotting = [] # Reset and populate with final state
    final_sample_to_cluster_map = {} # Create final map for reporting

    try:
        with open(cluster_file_path, "w") as f_cluster:
            # Sort clusters by number for consistent output
            sorted_cluster_names = sorted(current_clusters.keys(), key=lambda name: int(name.split('_')[-1]))

            for cluster_name in sorted_cluster_names:
                if cluster_name not in clusters_to_delete:
                    final_cluster_count += 1
                    samples_in_cluster = sorted(list(current_clusters[cluster_name]))
                    f_cluster.write(f"{cluster_name}\t{','.join(samples_in_cluster)}\n")
                    # Store info for plotting
                    cluster_info_for_plotting.append({
                        "name": cluster_name,
                        "samples": samples_in_cluster
                    })
                    # Update final map for reporting
                    for s in samples_in_cluster:
                        final_sample_to_cluster_map[s] = cluster_name

            if final_cluster_count == 0:
                 f_cluster.write("No linked clusters found.\n")
                 print("No final linked clusters to write.", file=sys.stderr)
            else:
                 print(f"Wrote {final_cluster_count} final cluster(s) to {cluster_file_path}", file=sys.stderr)
                 files_to_send.append(cluster_file_path)

    except IOError as e:
        print(f"Error writing final cluster file {cluster_file_path}: {e}", file=sys.stderr)

    # Note: final_sample_to_cluster_map now holds the definitive cluster assignment for each sample
    # The report generation logic later needs to be updated to use this map and the report_events dictionary.

    # --- Step 3b: Generate Cluster Plots (Looping after cluster file is written) ---
    print("\n--- Running Step 3b: Generating Cluster Plots ---", file=sys.stderr)
    if not cluster_info_for_plotting:
         print("No clusters found or generated, skipping plotting.", file=sys.stderr)
    else:
        for cluster_info in cluster_info_for_plotting:
            cluster_name = cluster_info["name"]
            sorted_samples = cluster_info["samples"]
            '''
            # --- Check if cluster needs processing ---
            if cluster_name not in modified_clusters_this_run:
                print(f"--- Skipping analysis for unchanged {cluster_name} ---", file=sys.stderr)
                # Still need to ensure the plot file is added to email list if it exists
                final_pdf_path = reports_dir / f"{cluster_name}_tsne.pdf"
                if final_pdf_path.is_file():
                     files_to_send.append(final_pdf_path)
                continue # Skip to the next cluster
            '''
            print(f"--- Analyzing and Generating Plot for {cluster_name} ---", file=sys.stderr)
            cluster_analysis_dir = None # Initialize for error handling before creation
            try:
                # Create a persistent directory for this cluster's analysis
                cluster_analysis_dir = args.base_dir / cluster_name
                cluster_analysis_dir.mkdir(parents=True, exist_ok=True)
                print(f"Using analysis directory for cluster: {cluster_analysis_dir}", file=sys.stderr)

                # 1. Combine FASTA files for cluster members
                cluster_combined_fasta = cluster_analysis_dir / f"{cluster_name}_combined.fasta" # Use specific name
                cluster_fasta_files = []
                with open(cluster_combined_fasta, "wb") as f_out_cat:
                    for sample in sorted_samples:
                        sample_fasta_gz = fasta_dir / f"{sample}.fasta.gz"
                        if sample_fasta_gz.is_file():
                            cluster_fasta_files.append(sample_fasta_gz)
                            try:
                                with gzip.open(sample_fasta_gz, "rb") as f_in:
                                    shutil.copyfileobj(f_in, f_out_cat)
                            except Exception as e: # Catch other MAFFT errors
                                print(f"Warning: Error reading or writing {sample_fasta_gz} for cluster: {e}", file=sys.stderr)
                        else:
                            print(f"Warning: FASTA file not found for cluster member {sample}: {sample_fasta_gz}", file=sys.stderr)

                if not cluster_fasta_files:
                    print(f"Error: No FASTA files found for any members of {cluster_name}. Skipping plot.", file=sys.stderr)
                    continue # Skip plot generation for this cluster

                # 2. Align combined FASTA using MAFFT (Docker)
                cluster_aligned_fasta = cluster_analysis_dir / f"{cluster_name}_aligned.fasta" # Use specific name
                mafft_docker_cmd = [
                    "docker", "run", "--rm",
                    "-v", f"{cluster_analysis_dir.resolve()}:/data", # Use resolved path for docker volume
                    "pegi3s/mafft",
                    "mafft", "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder",
                    f"/data/{cluster_combined_fasta.name}"
                ]
                print(f"Running MAFFT for {cluster_name}...", file=sys.stderr)
                mafft_ret = 1
                try:
                    # Log MAFFT command to stderr
                    print(f"Running MAFFT command for {cluster_name}: {' '.join(mafft_docker_cmd)}", file=sys.stderr)
                    with open(cluster_aligned_fasta, 'w') as f_out_aln:
                        proc_mafft = subprocess.run(mafft_docker_cmd, cwd=cluster_analysis_dir, text=True, capture_output=True, check=True)
                        f_out_aln.write(proc_mafft.stdout)
                        if proc_mafft.stderr:
                            # Print MAFFT stderr directly
                            print(f"MAFFT Stderr ({cluster_name}):\n{proc_mafft.stderr}", file=sys.stderr)
                    mafft_ret = 0
                except subprocess.CalledProcessError as e:
                    print(f"MAFFT for {cluster_name} failed (Exit Code {e.returncode}).", file=sys.stderr)
                    if e.stderr: print(f"MAFFT Error Output:\n{e.stderr}", file=sys.stderr) # Print error output
                except FileNotFoundError:
                    print(f"Error: Docker not found. Cannot run MAFFT for cluster plot.", file=sys.stderr)
                except Exception as e: # Catch other MAFFT errors
                    print(f"An unexpected error occurred running MAFFT for {cluster_name}: {e}", file=sys.stderr)

                if mafft_ret != 0 or not cluster_aligned_fasta.is_file() or cluster_aligned_fasta.stat().st_size == 0:
                    print(f"MAFFT failed or produced empty alignment for {cluster_name}. Skipping plot.", file=sys.stderr)
                    continue # Skip to next cluster

                # 3. Generate distance matrix using snp-dists
                cluster_pmatrix_file = cluster_analysis_dir / f"{cluster_name}_ident.pmatrix" # Use specific name
                cluster_dmatrix_file = cluster_analysis_dir / f"{cluster_name}_ident.dmatrix" # Use specific name
                # Use absolute paths in command string for robustness, use -m flag
                snp_dists_cmd_str = f"{args.snp_dists_path} -m {cluster_aligned_fasta.resolve()} > {cluster_pmatrix_file.resolve()}"
                snp_dists_cmd_str_dists = f"{args.snp_dists_path} -b {cluster_aligned_fasta.resolve()} > {cluster_dmatrix_file.resolve()}"
                print(f"Running snp-dists for {cluster_name}: {snp_dists_cmd_str}", file=sys.stderr)
                snp_dist_ret = 1
                try:
                    # Log snp-dists command to stderr
                    print(f"Running snp-dists command for {cluster_name}: {snp_dists_cmd_str}", file=sys.stderr)

                    proc_snp = subprocess.run(snp_dists_cmd_str, shell=True, cwd=cluster_analysis_dir, text=True, capture_output=True, check=False)
                    subprocess.run(snp_dists_cmd_str_dists, shell=True, cwd=cluster_analysis_dir, text=True, capture_output=False, check=False)
                    snp_dist_ret = proc_snp.returncode

                    if snp_dist_ret != 0:
                        print(f"snp-dists for {cluster_name} failed (Exit Code {snp_dist_ret}). Skipping plot.", file=sys.stderr)
                        if proc_snp.stderr:
                            print(f"Stderr:\n{proc_snp.stderr}", file=sys.stderr)
                        # Removed logging to file
                except Exception as e:
                    print(f"Error running snp-dists for {cluster_name}: {e}", file=sys.stderr)

                if snp_dist_ret != 0 or not cluster_pmatrix_file.is_file() or cluster_pmatrix_file.stat().st_size == 0:
                    print(f"snp-dists failed or produced empty matrix for {cluster_name}. Skipping plot.", file=sys.stderr)
                    continue # Skip to next cluster

                # 4. Run R script for t-SNE plot
                external_r_script_path = args.base_dir / "plot_tsne.R" # Assuming it's in base_dir
                if not external_r_script_path.is_file():
                     print(f"Error: R script plot_tsne.R not found at {external_r_script_path}. Cannot generate cluster plot.", file=sys.stderr)
                     continue # Skip plot for this cluster

                print(f"Running R script for {cluster_name} plot...", file=sys.stderr)
                # Log R command to stderr
                # Revert to original positional argument call
                r_cmd_list = [args.rscript_path, str(external_r_script_path.resolve()), str(cluster_pmatrix_file.resolve())]
                print(f"Running R command for {cluster_name}: {' '.join(r_cmd_list)}", file=sys.stderr)

                # Run R script directly using subprocess.run, not run_script helper
                r_ret = 1 # Default to error
                try:
                    # Execute Rscript, capture output
                    r_process = subprocess.run(r_cmd_list, cwd=cluster_analysis_dir, text=True, capture_output=True, check=False) # Run in cluster analysis dir
                    r_ret = r_process.returncode
                    if r_ret != 0:
                        print(f"Error: R script failed for {cluster_name} plot (exit code {r_ret}).", file=sys.stderr)
                        if r_process.stderr: print(f"R Stderr:\n{r_process.stderr}", file=sys.stderr)
                        if r_process.stdout: print(f"R Stdout:\n{r_process.stdout}", file=sys.stderr) # Print stdout too on error
                    else:
                         if r_process.stdout: print(f"R script stdout for {cluster_name}:\n{r_process.stdout}", file=sys.stderr) # Print stdout on success too
                except Exception as r_e:
                     print(f"Error executing R script for {cluster_name}: {r_e}", file=sys.stderr)
                     r_ret = 1 # Ensure failure state


                # 5. Move and rename plot (original logic)
                # Construct expected PDF name based on the pmatrix filename used as input for R
                temp_pdf_output_name = cluster_pmatrix_file.stem + "_tsne.pdf" # e.g., Cluster_1_ident_tsne.pdf
                temp_pdf_output_path = cluster_analysis_dir / temp_pdf_output_name # PDF is created in the analysis dir
                final_pdf_path = reports_dir / f"{cluster_name}_tsne.pdf"

                if r_ret == 0 and temp_pdf_output_path.is_file():
                    print(f"Cluster t-SNE plot generated: {temp_pdf_output_path.name}", file=sys.stderr)
                    try:
                        shutil.move(str(temp_pdf_output_path), str(final_pdf_path))
                        print(f"Moved cluster plot to: {final_pdf_path}", file=sys.stderr)
                        files_to_send.append(final_pdf_path) # Add plot to email list
                    except Exception as e:
                        print(f"Error moving cluster plot {temp_pdf_output_path} to {final_pdf_path}: {e}", file=sys.stderr)
                elif r_ret == 0 and not temp_pdf_output_path.is_file():
                     print(f"Error: R script finished successfully for {cluster_name} but output PDF not found at expected location: {temp_pdf_output_path}", file=sys.stderr)
                # Error message for r_ret != 0 already printed above
                else: # r_ret == 0 but file doesn't exist
                     print(f"Error: R script finished for {cluster_name} but output PDF not found: {temp_pdf_output_path}", file=sys.stderr)

            except Exception as cluster_e:
                 print(f"Error processing cluster {cluster_name} for plotting: {cluster_e}", file=sys.stderr)


    # --- Generate Run Report ---
    print("\n--- Generating Run Report ---", file=sys.stderr)
    run_report_file = reports_dir / f"Run_Report_{datetime.date.today()}.txt"
    try:
        # Create a map of sample -> final cluster details for easy lookup
        # final_sample_to_cluster_map is created in Step 3
        # report_events is created in Step 3

        with open(run_report_file, "w") as f_report:
            f_report.write(f"Pipeline Run Report - {datetime.datetime.now()}\n\n")

            processed_samples_in_report = set() # Keep track of samples added

            # Iterate through the samples based on the order they appear in run_report_lines
            sample_order = []
            temp_sample_lines = defaultdict(list)
            current_processing_sample = None
            for line in run_report_lines:
                 if line.startswith("--- Sample:"):
                      current_processing_sample = line.split(":")[-1].strip().split(" ")[0]
                      if current_processing_sample not in sample_order:
                           sample_order.append(current_processing_sample)
                      temp_sample_lines[current_processing_sample].append(line) # Add header
                      temp_sample_lines[current_processing_sample].append("") # Add blank line
                 elif current_processing_sample:
                      temp_sample_lines[current_processing_sample].append(line)

            # Now write report sample by sample
            for sample in sample_order:
                processed_samples_in_report.add(sample)
                # Write the collected status/link lines first
                f_report.write("\n".join(temp_sample_lines[sample]) + "\n")

                # Determine cluster status based on events and final mapping
                cluster_status_reported = False
                if sample in report_events:
                    # Prioritize reporting specific events like merges, additions, new formations
                    merge_event = next((e for e in report_events[sample] if "caused merge" in e), None)
                    added_event = next((e for e in report_events[sample] if "added to existing cluster" in e), None)
                    formed_event = next((e for e in report_events[sample] if "formed new cluster" in e), None)

                    final_cluster_name = final_sample_to_cluster_map.get(sample)
                    if not final_cluster_name: # Should not happen if event occurred, but safety check
                         f_report.write(f"\n{sample} involved in cluster event, but final cluster unknown.\n")
                         cluster_status_reported = True
                    else:
                         members = sorted(list(current_clusters.get(final_cluster_name, set())))
                         members_str = ", ".join(members)

                         if merge_event:
                              # Extract the original clusters if possible (complex, using simpler message for now)
                              # Example: "GP00092 involved in merge resulting in Cluster_1 -> GP00092, GP00093, ..."
                              f_report.write(f"\n{sample} involved in merge resulting in {final_cluster_name} -> {members_str}\n")
                              cluster_status_reported = True
                         elif added_event:
                              # Example: "GP00092 added to cluster Cluster_1 -> GP00092, GP00093, ..."
                              f_report.write(f"\n{sample} added to {final_cluster_name} -> {members_str}\n")
                              cluster_status_reported = True
                         elif formed_event:
                              # Example: "GP00092 formed new cluster Cluster_3 -> GP00092, GP00101"
                              f_report.write(f"\n{sample} formed new {final_cluster_name} -> {members_str}\n")
                              cluster_status_reported = True
                         # If only "linked within" events, treat as "already part of" below

                # If no specific event reported, check if it's just part of a cluster
                if not cluster_status_reported and sample in final_sample_to_cluster_map:
                    cluster_name = final_sample_to_cluster_map[sample]
                    members = sorted(list(current_clusters.get(cluster_name, set())))
                    members_str = ", ".join(members)
                    # Check if it was already in that cluster initially
                    if sample in initial_sample_to_cluster and initial_sample_to_cluster[sample] == cluster_name:
                         f_report.write(f"\n{sample} already part of {cluster_name} -> {members_str}\n")
                         cluster_status_reported = True
                    else:
                         # This case might occur if a sample was part of a cluster that got merged *into* this one,
                         # but the sample itself didn't trigger the merge link. Report as "part of".
                         f_report.write(f"\n{sample} part of final {cluster_name} -> {members_str}\n")
                         cluster_status_reported = True


                # If still no cluster status reported, check for failure/no links based on buffered lines
                if not cluster_status_reported:
                     buffered_lines = "\n".join(temp_sample_lines[sample])
                     if "Links found: No" in buffered_lines or "Status: Kmer Analysis SKIPPED" in buffered_lines or "FAILED" in buffered_lines:
                          f_report.write(f"\n{sample} has no links and is not part of any cluster.\n")
                          cluster_status_reported = True
                     # else: # Sample processed but somehow didn't end up in a cluster or have an event?
                     #      f_report.write(f"\n{sample} - final cluster status undetermined.\n")


                # Add separator
                f_report.write("\n=================================================\n\n")

            # Add entries for any processed samples missed by the loop (e.g., failed very early)
            all_processed_samples = set(sample_temp_dirs.keys()) # Get all samples attempted
            missed_samples = all_processed_samples - processed_samples_in_report
            for sample in sorted(list(missed_samples)):
                 f_report.write(f"--- Sample: {sample} ---\n\n")
                 # Try to find a status line for it in the original run_report_lines
                 status_line = next((line for line in run_report_lines if sample in line and ("FAILED" in line or "SKIPPED" in line)), None)
                 if status_line:
                      f_report.write(status_line.strip() + "\n")
                 else:
                      f_report.write(f"Status: Processing details incomplete (check individual report: {reports_dir / f'{sample}_report.out'}).\n")
                 f_report.write(f"\n{sample} has no links and is not part of any cluster.\n")
                 f_report.write("\n=================================================\n\n")

            # Remove final blank line added by loop
            f_report.seek(f_report.tell() - 2)
            f_report.truncate()

        # This block is now outside the 'with open' but inside the 'try'
        print(f"Run report saved to: {run_report_file}", file=sys.stderr)
        files_to_send.append(run_report_file)
    except IOError as e:
        print(f"Error writing run report file {run_report_file}: {e}", file=sys.stderr)


    # --- Step 4: Emailing (Optional) ---
    print("\n--- Pipeline Finished ---", file=sys.stderr)

    print("Reports generated in:", reports_dir, file=sys.stderr)


if __name__ == "__main__":
    # Need datetime for timestamp
    # import datetime # Already imported above
    main()