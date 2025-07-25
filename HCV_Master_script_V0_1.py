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
import csv # Added for parsing snp-dists output
import numpy as np
import io
import pandas as pd

def mst_tree(mst_output, snp_matrix_dest, args, degree): 
    print(f"  Running MST visualization...", file=sys.stderr)

    mst_cmd = [
        sys.executable,
        str(args.base_dir / "build_mst.py"),
        str(snp_matrix_dest),
        str(mst_output),
        "--min_degree", str(degree),
    ]
    print(f"  Running MST command: {' '.join(mst_cmd)}", file=sys.stderr)
    result = subprocess.run(mst_cmd, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        print(f"  Warning: MST generation returned code {result.returncode}", file=sys.stderr)
        if result.stderr:
            print(f"  MST stderr: {result.stderr}", file=sys.stderr)


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

def get_base_sample_name(full_seq_name: str) -> str | None:
    """Extracts the base sample name (e.g., 1700046311) from a full sequence name (e.g., _R_158_1700046311_1)."""
    parts = full_seq_name.split('_')

    # Handle _R_* prefixed names (e.g., _R_158_1700046311_1)
    if len(parts) >= 3 and parts[0] == '' and parts[1] == 'R' and parts[2].isdigit():
        # The sample ID is the part before the last underscore
        base_name = parts[-2] if len(parts) > 3 else full_seq_name
        if base_name.isdigit():  # Ensure it's a numeric sample ID
            return base_name

    # Original handling for other formats
    if len(parts) >= 2:
        if parts[0].isdigit() and parts[-1].isdigit() and len(parts) >= 3:
            return "_".join(parts[1:-1])
        elif len(parts) == 1:
            return parts[0]
        if len(parts) > 1 and parts[1].startswith("GP"):
            return parts[1]

    print(f"Warning: Could not reliably parse base name from '{full_seq_name}'. Using full name.", file=sys.stderr)
    return full_seq_name


def calculate_shared_haplotype_percentage(square_matrix_path: Path, s1_base: str, s2_base: str) -> tuple[float, int, int, float, int, int]:


    """
    Calculates the percentage of shared haplotypes (0 SNP distance) between two samples
    based on a snp-dists square matrix.
    Returns (percent_s1_to_s2, s1_total_genotypes, s1_shared_genotypes,
               percent_s2_to_s1, s2_total_genotypes, s2_shared_genotypes).
    """
    s1_sequences = defaultdict(list) # {base_sample_name: [full_seq_id1, full_seq_id2, ...]}
    s2_sequences = defaultdict(list)

    # Read the square matrix and populate sequence lists
   
    with open(square_matrix_path, 'r') as f_sq:
        reader = csv.reader(f_sq, delimiter='\t')
        header_cols = next(reader) # This is the first row, including the empty string at [0]
        
        # Map full sequence IDs to their column index in the data rows
        # The first element of header_cols is usually empty, so actual sequence IDs start from index 1
        # We will use this to map column index to sequence ID for populating matrix_data_dict
        
        matrix_data_dict = {} # {row_seq_id: {col_seq_id: distance}}
        for row in reader:
            row_seq_id = row[0] # First element of each row is the sequence ID
            matrix_data_dict[row_seq_id] = {}
            for col_idx, distance_str in enumerate(row[1:]): # Distances start from the second element
                col_seq_id = header_cols[col_idx + 1] # +1 because header_cols has an empty string at index 0
                matrix_data_dict[row_seq_id][col_seq_id] = distance_str

        # Populate s1_sequences and s2_sequences
        
        # Populate s1_sequences and s2_sequences from the header_cols
        for full_seq_id in header_cols[1:]: # Skip the first empty element
            base_name = get_base_sample_name(full_seq_id)
            if base_name == s1_base:
                s1_sequences[s1_base].append(full_seq_id)
            elif base_name == s2_base:
                s2_sequences[s2_base].append(full_seq_id)



    # Initialize percentages and counts
    percent_s1_to_s2 = 0.0
    s1_total_genotypes = 0
    s1_shared_genotypes = 0
    percent_s2_to_s1 = 0.0
    s2_total_genotypes = 0
    s2_shared_genotypes = 0
    # Calculate percent_s1_to_s2
    s1_total_genotypes = len(s1_sequences[s1_base])
    if s1_total_genotypes > 0:
        for s1_full_seq_id in s1_sequences[s1_base]:
            for s2_full_seq_id in s2_sequences[s2_base]:
                try:
                    distance = int(float(matrix_data_dict[s1_full_seq_id][s2_full_seq_id]))
                    if distance == 0:
                        s1_shared_genotypes += 1
                        break # Found at least one 0-SNP match for this s1 sequence
                except (ValueError, IndexError):
                    print(f"Warning: a Could not parse distance from matrix for ({s1_full_seq_id}, {s2_full_seq_id}). Value: '{matrix_data_dict.get(s1_full_seq_id, {}).get(s2_full_seq_id, 'N/A')}'", file=sys.stderr)
                    

    percent_s1_to_s2 = (s1_shared_genotypes / s1_total_genotypes) * 100 if s1_total_genotypes > 0 else 0.0

    # Calculate percent_s2_to_s1
    s2_total_genotypes = len(s2_sequences[s2_base])
    if s2_total_genotypes > 0:
        for s2_full_seq_id in s2_sequences[s2_base]:
            for s1_full_seq_id in s1_sequences[s1_base]:
                try:
                    distance = int(float(matrix_data_dict[s2_full_seq_id][s1_full_seq_id]))
                    if distance == 0:
                        s2_shared_genotypes += 1
                        break # Found at least one 0-SNP match for this s2 sequence
                except (ValueError, IndexError):
                    print(f"Warning: b Could not parse distance from matrix for ({s2_full_seq_id}, {s1_full_seq_id}). Value: '{matrix_data_dict.get(s2_full_seq_id, {}).get(s1_full_seq_id, 'N/A')}'", file=sys.stderr)
    return percent_s1_to_s2, s1_total_genotypes, s1_shared_genotypes, percent_s2_to_s1, s2_total_genotypes, s2_shared_genotypes


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
    ##parser.add_argument("--transmission_script", type=str, default="HCV_transmission_test_V0_2.py", help="Name of the transmission test script.")
    parser.add_argument("--r1_pattern", type=str, default="*_R1_001.fastq.gz", help="Glob pattern for R1 files.")
    parser.add_argument("--kmer_overlap_threshold", type=float, default=0.005, help="Threshold for k-mer overlap ratio to trigger transmission test.")
    parser.add_argument("--snp_dists_path", type=str, default="snp-dists", help="Path to the snp-dists executable.")
    parser.add_argument("--mafft_path", type=str, default="mafft", help="Path to the mafft executable.") # Added MAFFT path
    parser.add_argument("--rscript_path", type=str, default="Rscript", help="Path to the Rscript executable.")
    parser.add_argument("--hcv_kmer_threshold", type=float, default=0.1, help="Minimum fraction of matching kmers for a sample to be considered HCV")
    parser.add_argument("--snp_distance_threshold", type=float, default=9, help="SNP distance threshold for transmission confirmation.")
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


    reports_dir.mkdir(parents=True, exist_ok=True)

    combine_script_path = args.base_dir / args.combine_script
    distancer_script_path = args.base_dir / args.distancer_script

    if not combine_script_path.is_file():
         print(f"Error: Combine script not found: {combine_script_path}", file=sys.stderr)
         sys.exit(1)
    if not distancer_script_path.is_file():
         print(f"Error: Distancer script not found: {distancer_script_path}", file=sys.stderr)
         sys.exit(1)
    # Transmission script existence check moved to where it's needed

    processed_samples = []
    files_to_send = [] # For potential emailing
    run_report_lines = [] # Store lines for the final run report

    # --- Load Existing Clusters ---
    print("\n--- Loading Existing Clusters ---", file=sys.stderr)
    cluster_file_path = reports_dir / "clusters.txt" # Define the path to the cluster file
    existing_clusters_map, initial_sample_to_cluster, max_cluster_num, existing_links = load_existing_clusters(cluster_file_path)

    # Initialize current state for clustering
    current_clusters = existing_clusters_map.copy()  # {cluster_name: {samples}}
    current_sample_to_cluster = initial_sample_to_cluster.copy()  # {sample: cluster_name}
    next_cluster_id_counter = max_cluster_num + 1

    # Combine existing links with newly found links for downstream processing if needed, but clustering logic below uses only *new* links to modify state
    all_linked_pairs = existing_links.copy() # Start with links implied by existing clusters
    potential_kmer_links = set() # Track potential links found by kmer analysis in this run
    newly_linked_pairs_this_run = set() # Track SNP-confirmed links found *in this specific run* for clustering

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
    
    for sample_prefix in processed_samples:
        print(f"\nAnalyzing k-mers for sample: {sample_prefix}", file=sys.stderr)
        report_file = reports_dir / f"{sample_prefix}_report.out"
        kmer_file = kmer_dir / f"{sample_prefix}.kmers.gz"
        fasta_file = fasta_dir / f"{sample_prefix}.fasta.gz"
        temp_run_dir = sample_temp_dirs.get(sample_prefix)
        
        if not temp_run_dir:
             print(f"Error: Could not find temporary directory for sample {sample_prefix}. Skipping.", file=sys.stderr)
             run_report_lines.append(f"  Status: SKIPPED (Missing temp dir)")
             #sys.exit()
             continue
        
        try:
            with open(report_file, 'a') as f_out:
                f_out.write("\n--- HCV Reference Kmer Screening ---\n")
                try:
                    ref_dir = args.base_dir / "reference_genome"
                    is_hcv, virus_type, overlap1a, overlap1b = screen_sample_hcv(kmer_file, ref_dir, args.hcv_kmer_threshold)
                    f_out.write(f"Overlap with HCV 1a: {overlap1a*100:.2f}%\n")
                    f_out.write(f"Overlap with HCV 1b: {overlap1b*100:.2f}%\n")
                except Exception as e:
                    f_out.write(f"ERROR during HCV screening: {e}\n")
                    run_report_lines.append(f"{sample_prefix}: HCV screening failed with error: {e}")
                    #sys.exit()
                    continue
                f_out.write("--- End HCV Reference Kmer Screening ---\n")

                f_out.write("\n--- Kmer Distance Analysis ---\n")
                if not kmer_file.is_file():
                    f_out.write(f"ERROR: Kmer file not found: {kmer_file}\nSkipping distance analysis.\n")
                    print(f"Warning: Kmer file not found for {sample_prefix}, skipping distance analysis.", file=sys.stderr)
                    run_report_lines.append(f"  Status: Kmer Analysis SKIPPED (No kmer file)")
                    #sys.exit()
                    continue

                distancer_args = [str(args.base_dir), sample_prefix, str(temp_run_dir),  f"--overlap_threshold={args.kmer_overlap_threshold}"]
                print (f"Distancer arguments: {distancer_args}")
                print (f"Running kmer distancer script for {sample_prefix} with args: {distancer_args}", file=sys.stderr)
                return_code, stdout = run_script(distancer_script_path, distancer_args, capture_stdout=True)
                if stdout:
                    print(f"Kmer distancer script output for {sample_prefix}:\n{stdout}", file=sys.stderr)
                else:
                    print(f"No output from kmer distancer script for {sample_prefix}.", file=sys.stderr)
                if return_code != 0 or not stdout:
                    f_out.write(f"ERROR: Kmer distancer script failed with exit code {return_code}.\n")
                    print(f"ERROR: Kmer distancer failed for {sample_prefix}.", file=sys.stderr)
                    run_report_lines.append(f"  Status: Kmer Analysis FAILED (Code: {return_code})")
                    #sys.exit()
                    continue
                f_out.write(stdout)
                f_out.write("\n--- End Kmer Distance Analysis ---\n")

                # Parse kmer distancer output to identify linked samples
                distancer_output_file = temp_run_dir / f"{sample_prefix}_kmer_links.tsv"
                linked_samples = set()

                try:
                    with open(distancer_output_file, "r") as f_dist:
                        for line in f_dist:
                            if line.startswith("S1"):  # Skip the header line
                                continue
                            parts = line.strip().split()
                            if len(parts) == 6:  # Ensure the line has the expected format
                                sample1, sample2, overlap, total1, total2, ratio = parts
                                if float(ratio) >= args.kmer_overlap_threshold:
                                    # Add the newly found link to the set for this run's clustering update
                                    # Ensure consistent order by sorting the pair
                                    potential_link_pair = tuple(sorted((sample1, sample2)))
                                    potential_kmer_links.add(potential_link_pair) # Store as potential link
                                    # Add linked samples to the local set for reporting (optional, but was here)
                                    if sample1 == sample_prefix:
                                        linked_samples.add(sample2)
                                    elif sample2 == sample_prefix:
                                        linked_samples.add(sample1)
                except FileNotFoundError:
                    print(f"Warning: Kmer distancer output not found for {sample_prefix}. Skipping SNP confirmation for this sample.", file=sys.stderr)
                    run_report_lines.append(f"  Status: SKIPPED Kmer Link Parsing (File Not Found: {distancer_output_file})")
                    continue # Continue to the next sample
                except Exception as e:
                    print(f"Error parsing kmer distancer output for {sample_prefix}: {e}. Skipping SNP confirmation for this sample.", file=sys.stderr)
                    run_report_lines.append(f"  Status: FAILED Kmer Link Parsing (Error: {e})")
                    continue # Continue to the next sample
                # Report potential links found (moved reporting logic slightly)
                # Find partners linked *to this specific sample* in the potential links found so far
                partners_for_this_sample = {s for pair in potential_kmer_links for s in pair if s != sample_prefix and sample_prefix in pair}
                if partners_for_this_sample:
                     print(f"Potential kmer links involving {sample_prefix} found: {', '.join(sorted(list(partners_for_this_sample)))}", file=sys.stderr)
                # Note: The immediate clustering logic based on k-mers has been removed.
                # SNP confirmation will happen after all samples are processed.

        except Exception as e:
            print(f"Error during Step 2 (kmer/SNP check) for {sample_prefix}: {e}", file=sys.stderr)
            run_report_lines.append(f"  Status: FAILED during Kmer/SNP Check (Exception: {e})")
            if report_file.exists():
                 with open(report_file, 'a') as f_out:
                      f_out.write(f"\nFATAL ERROR during kmer/SNP check analysis: {e}\n")

        # Clean up will happen after sample analysis is complete
        finally:
            pass



    # --- Step 2.5: SNP Distance Confirmation via Group Alignment ---
    print("\n--- Running Step 2.5: SNP Distance Confirmation via Group Alignment ---", file=sys.stderr)
    newly_linked_pairs_this_run.clear() # Ensure it's empty before populating
    snp_confirmed_links_count = 0
    # Store generated distance matrices for potential reuse in plotting
    # Key: frozenset of samples in the group, Value: Path to distance matrix
    generated_distance_matrices = {}
    # Initialize sample_snp_distances before the try block
    sample_snp_distances = defaultdict(list)
    
    try:
        if not potential_kmer_links:
            print("No potential k-mer links found to confirm.", file=sys.stderr)
        else:
            # Find connected components (potential clusters) based on k-mer links
            potential_clusters = find_clusters(potential_kmer_links)
            print(f"Found {len(potential_clusters)} potential clusters based on k-mer links.", file=sys.stderr)

            for i, group_samples_set in enumerate(potential_clusters):
                group_samples = sorted(list(group_samples_set)) # Consistent ordering
                group_key = frozenset(group_samples) # Use frozenset as dict key
                group_name = f"potential_group_{i+1}"
                print(f"\nProcessing {group_name} with {len(group_samples)} samples: {', '.join(group_samples)}", file=sys.stderr)

                if len(group_samples) < 2:
                    print(f"  Skipping group {group_name}: needs at least 2 samples.", file=sys.stderr)
                    continue

                # Find a sample within this group that was processed in the current run to use its temp dir
                anchor_sample = None
                group_run_temp_dir = None
                for sample in group_samples:
                    if sample in sample_temp_dirs:
                        anchor_sample = sample
                        group_run_temp_dir = sample_temp_dirs[anchor_sample]
                        break # Use the first one found

                if not anchor_sample or not group_run_temp_dir or not group_run_temp_dir.is_dir():
                     print(f"  Error: Could not find an existing temp directory for any sample in group {group_name}. Samples processed this run: {list(sample_temp_dirs.keys())}. Skipping group.", file=sys.stderr)
                     continue
                print(f"  Using temp directory of sample '{anchor_sample}' ({group_run_temp_dir}) for group processing", file=sys.stderr)

                # Define file paths within the anchor sample's run temp dir
                combined_fasta_path = group_run_temp_dir / f"{group_name}_combined.fasta"
                aligned_fasta_path = group_run_temp_dir / f"{group_name}_aligned.fasta"
                distance_matrix_path = group_run_temp_dir / f"{group_name}_snp_distances.tsv"

                # 1. Combine FASTA files for the group
                print(f"  Combining FASTA files...", file=sys.stderr)
                valid_samples_in_group = []
                with open(combined_fasta_path, 'wb') as f_out:
                    for sample in group_samples:
                        fasta_gz = fasta_dir / f"{sample}.fasta.gz"
                        if fasta_gz.is_file():
                            with gzip.open(fasta_gz, 'rb') as f_in:
                                shutil.copyfileobj(f_in, f_out)
                            valid_samples_in_group.append(sample)
                        else:
                            print(f"  Warning: FASTA file not found for {sample}, excluding from group alignment.", file=sys.stderr)

                if len(valid_samples_in_group) < 2:
                    print(f"  Skipping group {group_name}: Fewer than 2 samples with FASTA files found.", file=sys.stderr)
                    continue

                # 2. Align the combined FASTA using MAFFT via Docker
                print(f"  Aligning combined FASTA using MAFFT (via Docker)...", file=sys.stderr)
                # MAFFT writes alignment to stdout, so we capture it
                # Mount the group's temp directory to /data in the container
                # Mount the sample's run temp directory, MAFFT reads/writes relative to that mount point
                mafft_docker_cmd = [
                    "docker", "run", "--rm",
                    "-v", f"{group_run_temp_dir.resolve()}:/data", # Mount the sample's run temp dir
                    "pegi3s/mafft", # Docker image
                    "mafft", "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder",
                    f"/data/{combined_fasta_path.name}" # Input file path inside container (/data/group_name_combined.fasta)
                ]
                print(f"  Running Docker MAFFT command: {' '.join(mafft_docker_cmd)}", file=sys.stderr)
                try:
                    # Run Docker command, capture stdout (alignment)
                    mafft_process_raw = subprocess.run(mafft_docker_cmd, check=True, text=True, capture_output=True)
                    mafft_process = mafft_process_raw.stdout.replace('_R_', '')
                    # Write captured stdout to the aligned fasta file
                    with open(aligned_fasta_path, 'w') as f_aligned_out:
                        f_aligned_out.write(mafft_process)

                    print(f"  MAFFT alignment complete: {aligned_fasta_path}", file=sys.stderr)
                except FileNotFoundError:
                     print(f"  Error: MAFFT executable not found at {args.mafft_path}. Skipping group {group_name}.", file=sys.stderr)
                     continue # Skip to next group
                except subprocess.CalledProcessError as e:
                     print(f"  Error running MAFFT for group {group_name}: {e}", file=sys.stderr)
                     if e.stderr: print(f"  MAFFT Stderr:\n{e.stderr}", file=sys.stderr)
                     continue # Skip to next group

                # 3. Run snp-dists on the aligned FASTA
                print(f"  Running snp-dists on aligned FASTA...", file=sys.stderr)
                # snp-dists needs the output file specified with -o
                snp_command = [args.snp_dists_path, "-m", "-a", str(aligned_fasta_path)] # Use -m for molten format, output goes to stdout
                try:
                    # Run snp-dists, capture stdout
                    snp_process = subprocess.run(snp_command, check=False, text=True, capture_output=True) # Use check=False to handle errors manually
                    if snp_process.returncode == 0:
                        # Write captured stdout to the distance matrix file
                        with open(distance_matrix_path, 'w') as f_dist_out:
                            f_dist_out.write(snp_process.stdout)
                        print(f"  snp-dists completed. Output written to: {distance_matrix_path}", file=sys.stderr)
                        # Store the path for potential reuse in plotting
                        generated_distance_matrices[group_key] = distance_matrix_path
                    else:
                        # Raise an error to be caught by the outer exception handler
                        raise subprocess.CalledProcessError(
                            returncode=snp_process.returncode,
                            cmd=snp_command,
                            stderr=snp_process.stderr,
                            stdout=snp_process.stdout
                        )
                    # Also generate the square matrix needed for plotting
                    square_matrix_path = group_run_temp_dir / f"{group_name}_snp_distances_square.tsv"
                    
                    snp_square_command = [args.snp_dists_path, "-b", "-a", str(aligned_fasta_path)] # Default square format to stdout
                    try:
                         snp_square_process = subprocess.run(snp_square_command, check=False, text=True, capture_output=True)
                         if snp_square_process.returncode == 0:
                              with open(square_matrix_path, 'w') as f_sq_out:
                                   f_sq_out.write(snp_square_process.stdout)
                              print(f"  snp-dists (square matrix for plotting) completed. Output: {square_matrix_path}", file=sys.stderr)
                              # Store the path to the square matrix for potential reuse in plotting
                              generated_distance_matrices[group_key] = square_matrix_path # Store path to SQUARE matrix


                    except Exception as plot_e:
                        # Catch errors during square matrix generation but don't stop link confirmation
                        print(f"  Warning: Error generating square snp-dists matrix for plotting: {plot_e}", file=sys.stderr)

                except FileNotFoundError:
                     print(f"  Error: snp-dists executable not found at {args.snp_dists_path}. Skipping group {group_name}.", file=sys.stderr)
                     continue # Skip to next group
                except subprocess.CalledProcessError as e:
                     print(f"  Error running snp-dists for group {group_name}: {e}", file=sys.stderr)
                     if e.stderr: print(f"  snp-dists Stderr:\n{e.stderr}", file=sys.stderr)
                     if e.stdout: print(f"  snp-dists Stdout:\n{e.stdout}", file=sys.stderr)
                     continue # Skip to next group

                # 4. Parse the 3-column distance output and confirm links
                print(f"  Parsing 3-column distance output to confirm links (Threshold <= {args.snp_distance_threshold})...", file=sys.stderr)
                
                # Initialize a dictionary to store SNP distances for reporting
                #sample_snp_distances = defaultdict(list)

                try:
                    with open(distance_matrix_path, 'r') as f_dist:
                        reader = csv.reader(f_dist, delimiter='\t')
                        
                        # New data structure to store all sequence-level distances
                        # Key: tuple of sorted base sample names, Value: list of distances
                        all_sequence_pair_distances = defaultdict(list)

                        for row in reader:
                            if len(row) == 3:
                                s1_full, s2_full, distance_str = row
                                s1_base = get_base_sample_name(s1_full)
                                s2_base = get_base_sample_name(s2_full)

                                if s1_base and s2_base and s1_base != s2_base: # Only consider inter-sample distances
                                    try:
                                        distance = int(float(distance_str))
                                        # Store the distance for later minimum calculation
                                        all_sequence_pair_distances[tuple(sorted((s1_base, s2_base)))].append(distance)
                                    except ValueError:
                                        print(f"    Warning: Could not parse distance value '{distance_str}' for pair ({s1_full}, {s2_full})", file=sys.stderr)
                                #else:
                                #     print(f"    Warning: Skipping malformed or intra-sample line in distance file: {row}", file=sys.stderr)
                                #     print(f"  Error: Malformed line in distance file: {row} x {distance_matrix_path}", file=sys.stderr)
                            elif row:
                                 print(f"    Warning: Skipping malformed line in distance file: {row}", file=sys.stderr)

                    # After parsing all rows, determine the minimum SNP distance for each sample pair
                    for (s1_base, s2_base), distances_list in all_sequence_pair_distances.items():
                        
                        if distances_list:
                            min_sample_distance = min(distances_list)                       
                            # Now apply the threshold and add to confirmed links
                            if min_sample_distance <= args.snp_distance_threshold:
                                confirmed_pair_base = tuple(sorted((s1_base, s2_base)))
                                if confirmed_pair_base not in newly_linked_pairs_this_run:
                                    print(f"    -> Link confirmed: {s1_base} <-> {s2_base} (SNP Distance: {min_sample_distance})", file=sys.stderr)
                                    newly_linked_pairs_this_run.add(confirmed_pair_base)
                                    snp_confirmed_links_count += 1
                                    all_linked_pairs.add(confirmed_pair_base)
                                    # Calculate shared haplotype percentage using the new method
                                    # The square_matrix_path is available from the group processing
                                    # Calculate shared haplotype percentage using the new method
                                    # The square_matrix_path is available from the group processing
                                    percent_s1_to_s2, s1_total_genotypes, s1_shared_genotypes, percent_s2_to_s1, s2_total_genotypes, s2_shared_genotypes =  calculate_shared_haplotype_percentage(square_matrix_path, s1_base, s2_base)
                                    
                                    # Store SNP distance and shared haplotype percentage for reporting
                                    # For s1_base, report percentage relative to s1_base's total genotypes
                                    sample_snp_distances[s1_base].append((s2_base, min_sample_distance, percent_s1_to_s2, s1_shared_genotypes, s1_total_genotypes))
                                    # For s2_base, report percentage relative to s2_base's total genotypes
                                    sample_snp_distances[s2_base].append((s1_base, min_sample_distance, percent_s2_to_s1, s2_shared_genotypes, s2_total_genotypes))

                                


                except FileNotFoundError:
                     print(f"  Error: Distance matrix file not found at {distance_matrix_path}. Cannot confirm links.", file=sys.stderr)
                except Exception as e:
                    print(f"  Error parsing distance file {distance_matrix_path} for group {group_name}: {e}", file=sys.stderr)

            print(f"\nSNP Confirmation finished. {snp_confirmed_links_count} links confirmed across all groups.", file=sys.stderr)

    except Exception as e:
         print(f"Error during SNP confirmation step: {e}", file=sys.stderr)
         # Decide if this should be fatal or just a warning
    # No finally block needed here; temp dirs are handled by the main loop's finally based on args.keep_tmp

    # --- Step 3: Stateful Clustering Update ---
    # This step now uses 'newly_linked_pairs_this_run' which only contains SNP-confirmed links from this run.
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
            '''
            snp_matrix_src = linked_matrix
            snp_matrix_dest = os.path.join(c1, f"{c1}_snp_distances_square.tsv")
            print (f"xXxXxXxXxX Source SNP distance matrix: {snp_matrix_src}", file=sys.stderr)
            aligned_fasta_path_dest = os.path.join(c1, f"{c1}_aligned.fasta")
            # Copy aligned FASTA to the new cluster directory   
            if os.path.exists(linked_fasta):
                shutil.copy2(linked_fasta, aligned_fasta_path_dest)
                print(f"Copied {linked_fasta} to {aligned_fasta_path_dest}", file=sys.stderr)
            # Copy SNP distance matrix and PDF to the new cluster directory
            if os.path.exists(snp_matrix_src):
                shutil.copy2(snp_matrix_src, snp_matrix_dest)
                print(f"Copied {snp_matrix_src} to {snp_matrix_dest}", file=sys.stderr)
            '''
        elif c2 and not c1 and c2 not in clusters_to_delete: # s2 in cluster, s1 is new
            print(f"Adding sample {s1} to cluster {c2} due to link with {s2}", file=sys.stderr)
            current_clusters[c2].add(s1)
            current_sample_to_cluster[s1] = c2
            modified_clusters_this_run.add(c2) # Mark cluster as modified
            report_events[s1].append(f"added to existing cluster {c2}")
            '''
            snp_matrix_src = linked_matrix
            snp_matrix_dest = os.path.join(c2, f"{c2}_snp_distances_square.tsv")
            print (f"xXxXxXxXxX Source SNP distance matrix: {snp_matrix_src}", file=sys.stderr)
            aligned_fasta_path_dest = os.path.join(c2, f"{c2}_aligned.fasta")
            # Copy aligned FASTA to the new cluster directory   
            if os.path.exists(linked_fasta):
                shutil.copy2(linked_fasta, aligned_fasta_path_dest)
                print(f"Copied {linked_fasta} to {aligned_fasta_path_dest}", file=sys.stderr)
            # Copy SNP distance matrix and PDF to the new cluster directory
            if os.path.exists(snp_matrix_src):
                shutil.copy2(snp_matrix_src, snp_matrix_dest)
                print(f"Copied {snp_matrix_src} to {snp_matrix_dest}", file=sys.stderr)
            '''
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
            # Create directory for the new cluster within Reports
            cluster_dir = new_cluster_name
            os.makedirs(cluster_dir, exist_ok=True)
            print(f"Created directory: {cluster_dir}", file=sys.stderr)
            
            # Get absolute path for cluster directory
            cluster_dir_abs = os.path.abspath(cluster_dir)
            
            #iterate over the samples in the new cluster and create the aligned FASTA and SNP distance matrix in the new cluster directory
            with open(os.path.join(cluster_dir, f"{new_cluster_name}_unaligned.fasta"), 'w') as f_out:
                for sample in new_cluster_set:
                    fasta_gz = fasta_dir / f"{sample}.fasta.gz"
                    if fasta_gz.is_file():
                        with gzip.open(fasta_gz, 'rt') as f_in:
                            shutil.copyfileobj(f_in, f_out)
                    else:
                        print(f"Warning: FASTA file not found for {sample}, skipping.", file=sys.stderr)
            # Align the combined FASTA using MAFFT via Docker
            print(f"Aligning combined FASTA for new cluster {new_cluster_name} using MAFFT (via Docker)...", file=sys.stderr)
            mafft_docker_cmd = [
                    "docker", "run", "--rm",
                    "-v", f"{cluster_dir_abs}:/data", # Use absolute path for mounting
                    "pegi3s/mafft",
                    "mafft", "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder",
                    f"/data/{new_cluster_name}_unaligned.fasta"
                ]
            print(f"Using Docker mount path: {cluster_dir_abs}:/data", file=sys.stderr)  # Debug line
            try:
                # Run Docker command, capture stdout (alignment)
                mafft_process_raw = subprocess.run(mafft_docker_cmd, check=True, text=True, capture_output=True)
                mafft_process = mafft_process_raw.stdout.replace('_R_', '')
                # Write captured stdout to the aligned fasta file
                with open(os.path.join(cluster_dir, f"{new_cluster_name}_aligned.fasta"), 'w') as f_aligned_out:
                    f_aligned_out.write(mafft_process)

                print(f"  MAFFT alignment complete: {os.path.join(cluster_dir, f'{new_cluster_name}_aligned.fasta')}", file=sys.stderr)
            except FileNotFoundError:
                 print(f"  Error: MAFFT executable not found at {args.mafft_path}. Skipping group {new_cluster_name}.", file=sys.stderr)
                 continue # Skip to next group
            except subprocess.CalledProcessError as e:
                 print(f"  Error running MAFFT for group {new_cluster_name}: {e}", file=sys.stderr)
                 if e.stderr: print(f"  MAFFT Stderr:\n{e.stderr}", file=sys.stderr)
                 continue # Skip to next group
            # create SNP distance matrix using snp-dists
            print(f"Running snp-dists on aligned FASTA for new cluster {new_cluster_name}...", file=sys.stderr)
            snp_command = [args.snp_dists_path, "-b", "-a", os.path.join(cluster_dir, f"{new_cluster_name}_aligned.fasta")]
            snp_process = subprocess.run(snp_command, check=False, text=True, capture_output=True) # Use check=False to handle errors manually
            # Parse output as 2D array (distance matrix)
            try:
                # Load SNP matrix as strings
                raw_matrix = np.genfromtxt(io.StringIO(snp_process.stdout), delimiter='\t', dtype=str)

                # Parse original sample names
                full_sample_names = raw_matrix[0, 1:]
                sample_base_names = []

                # Extract "biological" sample names between first and last underscore
                for name in full_sample_names:
                    parts = name.split('_')
                    if len(parts) >= 3:
                        sample_base_names.append('_'.join(parts[1:-1]))
                    else:
                        sample_base_names.append(name)  # Fallback: use full name if format unexpected

                full_sample_names = np.array(full_sample_names)
                sample_base_names = np.array(sample_base_names)

                # Convert matrix part to float
                dist_matrix = raw_matrix[1:, 1:].astype(float)

                # Fill diagonal to ignore self comparisons
                np.fill_diagonal(dist_matrix, np.inf)

                # Keep samples with distance <= threshold to a *different* base sample
                keep_mask = []
                for i, (sample_i, base_i) in enumerate(zip(full_sample_names, sample_base_names)):
                    row = dist_matrix[i]
                    base_diffs = (sample_base_names != base_i)
                    within_threshold = (row <= args.snp_distance_threshold)
                    should_keep = np.any(base_diffs & within_threshold)
                    keep_mask.append(should_keep)

                keep_mask = np.array(keep_mask)

                # Filter matrix and names
                filtered_matrix = dist_matrix[keep_mask][:, keep_mask]
                kept_sample_names = full_sample_names[keep_mask]

                # Reconstruct full matrix with headers and row labels
                output_lines = []
                output_lines.append('\t' + '\t'.join(kept_sample_names))
                for name, row in zip(kept_sample_names, filtered_matrix):
                    row_str = '\t'.join(str(int(val)) if np.isfinite(val) else '0' for val in row)
                    output_lines.append(f"{name}\t{row_str}")

                output_path = os.path.join(cluster_dir, f"{new_cluster_name}_snp_distances.tsv")
                with open(output_path, 'w') as f_dist_out:
                    f_dist_out.write('\n'.join(output_lines) + '\n')

                print(f"  Filtered SNP matrix written to: {output_path}", file=sys.stderr)

            except Exception as e:
                print(f"\n\n!!!! Error parsing or filtering SNP distance matrix for {new_cluster_name}: {e} !!!!\n\n", file=sys.stderr)
                continue


            # run MST visualization for the new cluster
            html_output = (os.path.join(cluster_dir, f"{new_cluster_name}_MST_full.html"))
            matrix_output = (os.path.join(cluster_dir, f"{new_cluster_name}_snp_distances.tsv"))
            print(f"Running MST visualization for new cluster {new_cluster_name}...", file=sys.stderr)
            # Call the mst_tree function to generate the MST visualization
            mst_tree(html_output, matrix_output, args, 0) # Use degree 0 for new clusters
            html_output = (os.path.join(cluster_dir, f"{new_cluster_name}_MST_overview.html"))
            mst_tree(html_output, matrix_output, args, 1) # Use degree 0 for new clusters

            next_cluster_id_counter += 1
    # --- Write the final cluster state to file ---
    cluster_file_path = reports_dir / "clusters.txt"
    try:
        with open(cluster_file_path, "w") as f_cluster:
            for cluster_name, samples in current_clusters.items():
                f_cluster.write(f"{cluster_name}\t{','.join(sorted(samples))}\n")
        print(f"Updated clusters written to {cluster_file_path}.", file=sys.stderr)
    except Exception as e:
        print(f"Error writing clusters to {cluster_file_path}: {e}", file=sys.stderr)

    # Note: final_sample_to_cluster_map now holds the definitive cluster assignment for each sample
    # The report generation logic later needs to be updated to use this map and the report_events dictionary.
    final_sample_to_cluster_map = current_sample_to_cluster.copy()

    # Prepare cluster information for report
    cluster_info_for_plotting = []
    # Use the final clusters after potential merges/deletions
    final_clusters_after_update = {k: v for k, v in current_clusters.items() if k not in clusters_to_delete}
    for cluster_name, cluster_samples in final_clusters_after_update.items():
         # Only plot clusters that were newly created or modified *structurally* this run
         # (i.e., added samples, merged, or newly formed)
         if cluster_name in modified_clusters_this_run:
              cluster_info_for_plotting.append({
                   "name": cluster_name,
                   "samples": cluster_samples # Use the final set of samples
              })

    # --- Generate Run Report ---
    print("\n--- Generating Run Report ---", file=sys.stderr)
    run_report_file = reports_dir / f"Run_Report_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.txt"

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
            # Add SNP distance details to the report
            if sample in sample_snp_distances and sample_snp_distances[sample]:
                # Sort distances for consistent reporting
                sorted_distances = sorted(sample_snp_distances[sample], key=lambda item: item[1])
                
                # Extract total_count_for_sample once from the first linked pair (all should be the same for 'sample')
                # The tuple is (linked_sample, distance, shared_haplotype_percent, shared_count, total_count_for_sample)
                if sorted_distances: # Ensure there's at least one linked sample
                    _, _, _, _, total_haplotypes_in_sample = sorted_distances[0]
                    f_report.write(f"\nTotal haplotypes in {sample}: {total_haplotypes_in_sample}\n")

                f_report.write(f"\n  Minimum SNP Distances to Other Samples (<= {args.snp_distance_threshold}):\n")
                for linked_sample, distance, shared_haplotype_percent, shared_count, _ in sorted_distances: # Ignore total_count_for_sample here
                    # Add the shared haplotype percentage and raw count to the report line
                    shared_haplotype_percent = (shared_count / total_haplotypes_in_sample * 100) if total_haplotypes_in_sample > 0 else 0
                    f_report.write(f"    -> Link confirmed: {sample} <-> {linked_sample} (SNP Distance: {distance}), {shared_haplotype_percent:.2f}% shared haplotypes ({shared_count} haplotypes)\n")
                    
                f_report.write("\n") # Add a blank line after the distances

            # Determine cluster status based on events and final mapping
            cluster_status_reported = False
            if sample in report_events:
                merge_event = next((e for e in report_events[sample] if "caused merge" in e), None)
                added_event = next((e for e in report_events[sample] if "added to existing cluster" in e), None)
                formed_event = next((e for e in report_events[sample] if "formed new cluster" in e), None)

                final_cluster_name = final_sample_to_cluster_map.get(sample)
                if final_cluster_name:
                        # Resolve final cluster name in case of chained merges
                        while final_cluster_name in merged_cluster_log:
                            final_cluster_name = merged_cluster_log[final_cluster_name]
                        members = sorted(list(final_clusters_after_update.get(final_cluster_name, set())))
                        members_str = ", ".join(members)
                        if merge_event:
                            f_report.write(f"\n{sample} involved in merge resulting in {final_cluster_name} -> {members_str}\n")
                            cluster_status_reported = True
                        elif added_event:
                            f_report.write(f"\n{sample} added to {final_cluster_name} -> {members_str}\n")
                            cluster_status_reported = True
                        elif formed_event:
                            f_report.write(f"\n{sample} formed new {final_cluster_name} -> {members_str}\n")
                            cluster_status_reported = True
                else: # Should not happen if event occurred, but safety check
                        f_report.write(f"\n{sample} involved in cluster event, but final cluster unknown.\n")
                        cluster_status_reported = True

            # If no specific event reported, check if it's just part of a cluster
            if not cluster_status_reported and sample in final_sample_to_cluster_map:
                cluster_name = final_sample_to_cluster_map[sample]
                # Resolve final cluster name in case of chained merges
                while cluster_name in merged_cluster_log:
                        cluster_name = merged_cluster_log[cluster_name]
                members = sorted(list(final_clusters_after_update.get(cluster_name, set())))
                members_str = ", ".join(members)
                # Check if it was already in that cluster initially
                initial_cluster = initial_sample_to_cluster.get(sample)
                while initial_cluster in merged_cluster_log: # Resolve initial cluster too
                        initial_cluster = merged_cluster_log[initial_cluster]

                if initial_cluster == cluster_name:
                        f_report.write(f"\n{sample} already part of {cluster_name} -> {members_str}\n")
                else: # Part of final cluster, but wasn't initially or added/merged explicitly
                        f_report.write(f"\n{sample} part of final {cluster_name} -> {members_str}\n")
                cluster_status_reported = True

            # If still no cluster status reported, assume no links/not clustered
            if not cluster_status_reported:
                    f_report.write(f"\n{sample} has no confirmed links and is not part of any cluster.\n")

            f_report.write("\n" + "="*50 + "\n\n") # Separator

        # Add entries for any processed samples missed by the loop (e.g., failed very early)
        all_processed_samples = set(sample_temp_dirs.keys()) # Get all samples attempted
        missed_samples = all_processed_samples - processed_samples_in_report
        if missed_samples:
                f_report.write("Samples with Incomplete Processing:\n")
                f_report.write("-" * 30 + "\n")
                for sample in sorted(list(missed_samples)):
                    f_report.write(f"--- Sample: {sample} ---\n")
                    # Try to find a status line for it in the original run_report_lines
                    status_line = next((line for line in run_report_lines if sample in line and ("FAILED" in line or "SKIPPED" in line)), None)
                    if status_line:
                        f_report.write(status_line.strip() + "\n")
                    else:
                        f_report.write(f"Status: Processing details incomplete (check individual report: {reports_dir / f'{sample}_report.out'}).\n")
                    f_report.write(f"{sample} has no confirmed links and is not part of any cluster.\n")
                    f_report.write("\n" + "="*50 + "\n\n")


    print(f"Run report saved to: {run_report_file}", file=sys.stderr)
    files_to_send.append(run_report_file)



    # --- Append SNP Distance Details to Individual Reports ---
    print("\n--- Appending SNP Distance Details to Individual Reports ---", file=sys.stderr)
    if 'processed_samples' in locals() and 'sample_snp_distances' in locals():
        for sample in processed_samples:
            print (f"  Appending SNP distance details for {sample}...", file=sys.stderr)
            report_file = reports_dir / f"{sample}_report.out"

            if report_file.is_file() and sample in sample_snp_distances and sample_snp_distances[sample]:
                try:
                    with open(report_file, 'a') as f_out:
                        # Sort distances for consistent reporting
                        sorted_distances = sorted(sample_snp_distances[sample], key=lambda item: item[1])

                        # Extract total_count_for_sample once from the first linked pair
                        if sorted_distances:
                            _, _, _, _, total_haplotypes_in_sample = sorted_distances[0]
                            f_out.write(f"\nTotal haplotypes in {sample}: {total_haplotypes_in_sample}\n")

                        f_out.write(f"\nMinimum SNP Distances to Other Samples (<= {args.snp_distance_threshold}):\n")
                        for linked_sample, distance, shared_haplotype_percent, shared_count, _ in sorted_distances: # Ignore total_count_for_sample here
                            shared_haplotype_percent = (shared_count / total_haplotypes_in_sample * 100) if total_haplotypes_in_sample > 0 else 0
                            f_out.write(f"  -> Link confirmed: {sample} <-> {linked_sample} (SNP Distance: {distance}), {shared_haplotype_percent:.2f}% shared haplotypes ({shared_count} haplotypes)\n")
                        f_out.write("\n") # Add a blank line after the distances

                        # Determine cluster status based on events and final mapping for individual report
                        cluster_status_reported = False
                        if sample in report_events:
                            merge_event = next((e for e in report_events[sample] if "caused merge" in e), None)
                            added_event = next((e for e in report_events[sample] if "added to existing cluster" in e), None)
                            formed_event = next((e for e in report_events[sample] if "formed new cluster" in e), None)

                            final_cluster_name = final_sample_to_cluster_map.get(sample)
                            if final_cluster_name:
                                    # Resolve final cluster name in case of chained merges
                                    while final_cluster_name in merged_cluster_log:
                                        final_cluster_name = merged_cluster_log[final_cluster_name]
                                    members = sorted(list(final_clusters_after_update.get(final_cluster_name, set())))
                                    members_str = ", ".join(members)
                                    if merge_event:
                                        f_out.write(f"\n{sample} involved in merge resulting in {final_cluster_name} -> {members_str}\n")
                                        cluster_status_reported = True
                                    elif added_event:
                                        f_out.write(f"\n{sample} added to {final_cluster_name} -> {members_str}\n")
                                        cluster_status_reported = True
                                    elif formed_event:
                                        f_out.write(f"\n{sample} formed new {final_cluster_name} -> {members_str}\n")
                                        cluster_status_reported = True
                            else: # Should not happen if event occurred, but safety check
                                    f_out.write(f"\n{sample} involved in cluster event, but final cluster unknown.\n")
                                    cluster_status_reported = True

                        # If no specific event reported, check if it's just part of a cluster
                        if not cluster_status_reported and sample in final_sample_to_cluster_map:
                            cluster_name = final_sample_to_cluster_map[sample]
                            # Resolve final cluster name in case of chained merges
                            while cluster_name in merged_cluster_log:
                                    cluster_name = merged_cluster_log[cluster_name]
                            members = sorted(list(final_clusters_after_update.get(cluster_name, set())))
                            members_str = ", ".join(members)
                            # Check if it was already in that cluster initially
                            initial_cluster = initial_sample_to_cluster.get(sample)
                            while initial_cluster in merged_cluster_log: # Resolve initial cluster too
                                    initial_cluster = merged_cluster_log[initial_cluster]

                            if initial_cluster == cluster_name:
                                    f_out.write(f"\n{sample} already part of {cluster_name} -> {members_str}\n")
                            else: # Part of final cluster, but wasn't initially or added/merged explicitly
                                    f_out.write(f"\n{sample} part of final {cluster_name} -> {members_str}\n")
                            cluster_status_reported = True

                        # If still no cluster status reported, assume no links/not clustered
                        if not cluster_status_reported:
                                f_out.write(f"\n{sample} has no confirmed links and is not part of any cluster.\n")
                        f_out.write("\n" + "="*50 + "\n\n") # Separator for individual report
                    print(f"  Appended SNP distance details to {report_file}", file=sys.stderr)
                except IOError as e:
                    print(f"  Error appending SNP distance details to {report_file}: {e}", file=sys.stderr)
                except Exception as e:
                    print(f"  Unexpected error appending SNP distance details for {sample}: {e}", file=sys.stderr)
            elif not report_file.is_file():
                    print(f"  Warning: Individual report file not found for {sample} at {report_file}. Cannot append SNP details.", file=sys.stderr)
            # else: sample_snp_distances[sample] is empty, no details to append

    else:
        print("  Warning: processed_samples or sample_snp_distances not found. Skipping appending SNP details to individual reports.", file=sys.stderr)


    # Clean up temp directories after all processing is complete
    for sample_prefix, temp_run_dir in sample_temp_dirs.items():
        if temp_run_dir and temp_run_dir.exists() and not args.keep_tmp:
            print(f"Removing temporary directory for sample {sample_prefix}: {temp_run_dir}", file=sys.stderr)
            shutil.rmtree(temp_run_dir, ignore_errors=True)
        elif temp_run_dir and temp_run_dir.exists() and args.keep_tmp:
            print(f"Keeping temporary directory for sample {sample_prefix}: {temp_run_dir}", file=sys.stderr)

    # --- Step 4: Emailing (Optional) ---
    print("\n--- Pipeline Finished ---", file=sys.stderr)

    print("Reports generated in:", reports_dir, file=sys.stderr)


if __name__ == "__main__":
    # Need datetime for timestamp
    # import datetime # Already imported above
    main()
