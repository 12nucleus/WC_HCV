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

# New helper functions for HCV reference screening
def load_reference_kmers(ref_file: Path) -> set:
    kmers = set()
    with open(ref_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                kmers.add(parts[0])
    return kmers

def load_sample_kmers(sample_kmer_file: Path) -> set:
    kmers = set()
    with gzip.open(sample_kmer_file, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                kmers.add(parts[0])
    return kmers

def screen_sample_hcv(sample_kmer_file: Path, ref_dir: Path, threshold: float):
    # Load reference kmers
    ref_file1 = ref_dir / "1a_kmer_ref_counts.out"
    ref_file2 = ref_dir / "1b_kmer_ref_counts.out"
    if not (ref_file1.is_file() and ref_file2.is_file()):
        raise FileNotFoundError(f"Reference kmer files not found in {ref_dir}")
    sample_kmers = load_sample_kmers(sample_kmer_file)
    ref1a_kmers = load_reference_kmers(ref_file1)
    ref1b_kmers = load_reference_kmers(ref_file2)
    overlap1a = len(sample_kmers & ref1a_kmers) / len(ref1a_kmers) if ref1a_kmers else 0
    overlap1b = len(sample_kmers & ref1b_kmers) / len(ref1b_kmers) if ref1b_kmers else 0
    is_hcv = (overlap1a >= threshold) or (overlap1b >= threshold)
    virus_type = "1a" if overlap1a >= overlap1b else "1b"
    return is_hcv, virus_type, overlap1a, overlap1b

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
    parser.add_argument("--hcv_kmer_threshold", type=float, default=0.1,
                        help="Minimum fraction of matching kmers for a sample to be considered HCV")
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
                f_out.write("\n--- HCV Reference Kmer Screening ---\n")
                try:
                    ref_dir = args.base_dir / "reference_genome"
                    is_hcv, virus_type, overlap1a, overlap1b = screen_sample_hcv(kmer_file, ref_dir, args.hcv_kmer_threshold)
                    f_out.write(f"Overlap with HCV 1a: {overlap1a*100:.2f}%\n")
                    f_out.write(f"Overlap with HCV 1b: {overlap1b*100:.2f}%\n")
                    #if not is_hcv:
                    #    f_out.write("Sample does not appear to be HCV based on reference kmer screening. Skipping further analysis.\n")
                    #    run_report_lines.append(f"{sample_prefix}: Not HCV based on kmer screening")
                    #    continue  # Skip to next sample
                    #else:
                    #    f_out.write(f"Sample determined to be HCV type {virus_type} based on reference kmer screening.\n")
                except Exception as e:
                    f_out.write(f"ERROR during HCV screening: {e}\n")
                    run_report_lines.append(f"{sample_prefix}: HCV screening failed with error: {e}")
                    continue
                f_out.write("--- End HCV Reference Kmer Screening ---\n")

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

                f_out.write(f"\n--- Combined SNP Distance Check (Threshold <= 9) ---\n")
                if not transmission_candidates:
                    f_out.write(f"No potential transmission candidates found based on k-mer overlap >= {args.kmer_overlap_threshold}.\n")
                    run_report_lines.append(f"  Links found: No (No k-mer candidates)")
                else:
                    # Build a list of FASTA files: current sample plus all subject candidates
                    fasta_list = [fasta_file]
                    missing_files = False
                    for subject_prefix in transmission_candidates:
                        subject_fasta_file = fasta_dir / f"{subject_prefix}.fasta.gz"
                        if subject_fasta_file.is_file():
                            fasta_list.append(subject_fasta_file)
                        else:
                            f_out.write(f"WARNING: Subject FASTA file not found: {subject_fasta_file}. Skipping this candidate.\n")
                    
                    if len(fasta_list) < 2:
                        f_out.write("Not enough FASTA files available for a combined SNP check.\n")
                        run_report_lines.append("  Links found: SKIPPED (Not enough samples for SNP check)")
                    else:
                        # Combine all chosen FASTA files into one
                        combined_fasta = temp_run_dir / f"combined_{sample_prefix}_all.fasta"
                        with open(combined_fasta, "wb") as f_cat:
                            for f_gz in fasta_list:
                                with gzip.open(f_gz, "rb") as f_in:
                                    shutil.copyfileobj(f_in, f_cat)
                        f_out.write(f"Combined FASTA created: {combined_fasta}\n")
                        
                        # Run MAFFT (via docker) on the combined FASTA
                        combined_aligned_fasta = temp_run_dir / f"combined_{sample_prefix}_all_aligned.fasta"
                        mafft_docker_cmd = [
                            "docker", "run", "--rm", "-v", f"{temp_run_dir.resolve()}:/data",
                            "pegi3s/mafft", "mafft",
                            "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder",
                            f"/data/{combined_fasta.name}"
                        ]
                        try:
                            with open(combined_aligned_fasta, 'w') as f_out_aln:
                                proc_mafft = subprocess.run(
                                    mafft_docker_cmd,
                                    cwd=temp_run_dir,
                                    text=True,
                                    capture_output=True,
                                    check=True
                                )
                                f_out_aln.write(proc_mafft.stdout)
                            f_out.write(f"MAFFT alignment completed: {combined_aligned_fasta}\n")
                        except Exception as e:
                            f_out.write(f"ERROR: MAFFT alignment failed for combined samples: {e}\n")
                            run_report_lines.append(f"  Links found: FAILED (MAFFT alignment error)")
                            continue

                        # Run snp-dists on the aligned FASTA
                        combined_snp_output = temp_run_dir / f"snp_distances_{sample_prefix}_combined.tsv"
                        snp_dists_cmd_str = f"{args.snp_dists_path} -m {combined_aligned_fasta.name} > {combined_snp_output.name}"
                        try:
                            proc_snp = subprocess.run(
                                snp_dists_cmd_str,
                                shell=True,
                                cwd=temp_run_dir,
                                text=True,
                                capture_output=True,
                                check=True
                            )
                            f_out.write(f"snp-dists executed successfully. Output file: {combined_snp_output}\n")
                        except Exception as e:
                            f_out.write(f"ERROR: snp-dists failed for combined samples: {e}\n")
                            run_report_lines.append("  Links found: FAILED (snp-dists error)")
                            continue

                        # Parse the SNP distance output (three column format: sample1 sample2 distance)
                        links = []
                        try:
                            with open(combined_snp_output, "r") as f_snp:
                                for line in f_snp:
                                    parts = line.strip().split()
                                    if len(parts) != 3:
                                        continue  # Skip malformed lines
                                    sample1, sample2, snp_val = parts
                                    try:
                                        distance = int(snp_val)
                                        if sample1 != sample2 and distance <= 9:
                                            links.append((sample1, sample2, distance))
                                    except ValueError:
                                        f_out.write(f"Warning: Could not parse SNP distance value: {snp_val}\n")
                        except Exception as e:
                            f_out.write(f"ERROR parsing snp-dists output: {e}\n")
                        
                        # Write SNP links to the sample report and add summary to run_report_lines
                        if links:
                            f_out.write(f"Found {len(links)} SNP link(s) (≤9 SNP differences):\n")

                            # Extract sample names from sequence IDs (e.g. 1_GP00098_8049 -> GP00098)
                            def extract_sample_name(seq_id):
                                # Assumes format: <number>_<sample>_<number>
                                parts = seq_id.split('_')
                                if len(parts) >= 2:
                                    return parts[1]
                                return seq_id

                            linked_samples = set()
                            for sample1, sample2, distance in links:
                                # If either sample1 or sample2 is the analyzed sample, add the other sample's name
                                s1 = extract_sample_name(sample1)
                                s2 = extract_sample_name(sample2)
                                if s1 == sample_prefix and s2 != sample_prefix:
                                    linked_samples.add(s2)
                                elif s2 == sample_prefix and s1 != sample_prefix:
                                    linked_samples.add(s1)

                            if linked_samples:
                                f_out.write(f"Samples with ≤9 SNPs from {sample_prefix}: {', '.join(sorted(linked_samples))}\n")
                                run_report_lines.append(f"{sample_prefix} clusters with: {', '.join(sorted(linked_samples))}")
                            else:
                                f_out.write(f"No samples with ≤9 SNPs from {sample_prefix}.\n")
                                run_report_lines.append(f"{sample_prefix} has no SNP links ≤9")
                        else:
                            f_out.write("No SNP links (≤9 differences) found among combined samples.\n")
                            run_report_lines.append("No SNP links (≤9 differences)")

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
                 #f_cluster.write("No linked clusters found.\n")
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
            
            # Only process cluster if it was modified
            if cluster_name not in modified_clusters_this_run:
                print(f"--- Skipping analysis for unchanged {cluster_name} ---", file=sys.stderr)
                final_pdf_path = reports_dir / f"{cluster_name}_tsne.pdf"
                if final_pdf_path.is_file():
                     files_to_send.append(final_pdf_path)
                continue  # Skip to the next cluster
            
            print(f"--- Analyzing and Generating Plot for {cluster_name} ---", file=sys.stderr)
            cluster_analysis_dir = None  # Initialize for error handling before creation
            try:
                # Create a persistent directory for this cluster's analysis
                cluster_analysis_dir = args.base_dir / cluster_name
                cluster_analysis_dir.mkdir(parents=True, exist_ok=True)
                print(f"Using analysis directory for cluster: {cluster_analysis_dir}", file=sys.stderr)

                # 1. Combine FASTA files for cluster members
                cluster_combined_fasta = cluster_analysis_dir / f"{cluster_name}_combined.fasta"
                cluster_fasta_files = []
                with open(cluster_combined_fasta, "wb") as f_out_cat:
                    for sample in sorted_samples:
                        sample_fasta_gz = fasta_dir / f"{sample}.fasta.gz"
                        if sample_fasta_gz.is_file():
                            cluster_fasta_files.append(sample_fasta_gz)
                            try:
                                with gzip.open(sample_fasta_gz, "rb") as f_in:
                                    shutil.copyfileobj(f_in, f_out_cat)
                            except Exception as e:
                                print(f"Warning: Error reading/writing {sample_fasta_gz} for {cluster_name}: {e}", file=sys.stderr)
                        else:
                            print(f"Warning: FASTA file not found for member {sample}: {sample_fasta_gz}", file=sys.stderr)

                if not cluster_fasta_files:
                    print(f"Error: No FASTA files for any members of {cluster_name}. Skipping plot.", file=sys.stderr)
                    continue

                # 2. Align combined FASTA using MAFFT (via Docker)
                cluster_aligned_fasta = cluster_analysis_dir / f"{cluster_name}_aligned.fasta"
                mafft_docker_cmd = [
                    "docker", "run", "--rm",
                    "-v", f"{cluster_analysis_dir.resolve()}:/data",
                    "pegi3s/mafft",
                    "mafft", "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder",
                    f"/data/{cluster_combined_fasta.name}"
                ]
                print(f"Running MAFFT for {cluster_name}...", file=sys.stderr)
                try:
                    with open(cluster_aligned_fasta, 'w') as f_out_aln:
                        proc_mafft = subprocess.run(mafft_docker_cmd, cwd=cluster_analysis_dir, text=True,
                                                      capture_output=True, check=True)
                        f_out_aln.write(proc_mafft.stdout)
                        if proc_mafft.stderr:
                            print(f"MAFFT stderr for {cluster_name}:\n{proc_mafft.stderr}", file=sys.stderr)
                except Exception as e:
                    print(f"ERROR: MAFFT failed for {cluster_name}: {e}", file=sys.stderr)
                    continue

                if not cluster_aligned_fasta.is_file() or cluster_aligned_fasta.stat().st_size == 0:
                    print(f"MAFFT produced empty alignment for {cluster_name}. Skipping plot.", file=sys.stderr)
                    continue

                # 3. Generate distance matrix using snp-dists
                cluster_pmatrix_file = cluster_analysis_dir / f"{cluster_name}_ident.pmatrix"
                snp_dists_cmd_str = f"{args.snp_dists_path} -m {cluster_aligned_fasta.resolve()} > {cluster_pmatrix_file.resolve()}"
                print(f"Running snp-dists for {cluster_name}: {snp_dists_cmd_str}", file=sys.stderr)
                try:
                    proc_snp = subprocess.run(snp_dists_cmd_str, shell=True, cwd=cluster_analysis_dir,
                                              text=True, capture_output=True, check=True)
                except Exception as e:
                    print(f"ERROR: snp-dists failed for {cluster_name}: {e}", file=sys.stderr)
                    continue

                if not cluster_pmatrix_file.is_file() or cluster_pmatrix_file.stat().st_size == 0:
                    print(f"snp-dists produced empty matrix for {cluster_name}. Skipping plot.", file=sys.stderr)
                    continue

                # 4. Run R script to generate t-SNE plot
                external_r_script_path = args.base_dir / "plot_tsne.R"
                if not external_r_script_path.is_file():
                     print(f"Error: R script plot_tsne.R not found at {external_r_script_path}. Skipping plot.", file=sys.stderr)
                     continue

                r_cmd_list = [args.rscript_path, str(external_r_script_path.resolve()), str(cluster_pmatrix_file.resolve())]
                print(f"Running R script for {cluster_name}: {' '.join(r_cmd_list)}", file=sys.stderr)
                try:
                    r_process = subprocess.run(r_cmd_list, cwd=cluster_analysis_dir, text=True,
                                               capture_output=True, check=False)
                    if r_process.returncode != 0:
                        print(f"ERROR: R script for {cluster_name} failed (exit code {r_process.returncode}).", file=sys.stderr)
                        continue
                except Exception as r_e:
                     print(f"Error executing R script for {cluster_name}: {r_e}", file=sys.stderr)
                     continue

                # 5. Move and rename the plot PDF
                temp_pdf_output_path = cluster_analysis_dir / (cluster_pmatrix_file.stem + "_tsne.pdf")
                final_pdf_path = reports_dir / f"{cluster_name}_tsne.pdf"
                if temp_pdf_output_path.is_file():
                    try:
                        shutil.move(str(temp_pdf_output_path), str(final_pdf_path))
                        print(f"Moved TSNE plot for {cluster_name} to: {final_pdf_path}", file=sys.stderr)
                        files_to_send.append(final_pdf_path)
                    except Exception as e:
                        print(f"Error moving plot for {cluster_name}: {e}", file=sys.stderr)
                else:
                    print(f"Error: Expected TSNE plot not found for {cluster_name} at {temp_pdf_output_path}", file=sys.stderr)

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
