# HCV_transmission_test_V0_3_snp_dists.py
# Uses snp-dists upfront for all pairwise distance calculations.
import argparse
import subprocess
import gzip
import shutil
import tempfile
from pathlib import Path
import sys
import os
import math
import numpy as np # Keep numpy for mean/std dev/Z-score calculations
import pandas as pd # Use pandas for easier matrix parsing
from Bio import SeqIO # Still needed for header manipulation

# Removed the calculate_pairwise_distances_optimized function

def run_command(cmd_list, cwd=None, log_file=None, shell=False):
    """Runs an external command, logs output, checks return code."""
    cmd_str = ' '.join(cmd_list) if not shell else cmd_list # Use list for non-shell, string for shell
    print(f"Running command: {cmd_str}", file=sys.stderr)
    stdout_pipe = subprocess.PIPE
    stderr_pipe = subprocess.PIPE
    if log_file:
        log_path = Path(log_file)
        try:
            stdout_handle = open(log_path, 'a')
            stderr_handle = stdout_handle
            stdout_handle.write(f"\n--- Running: {cmd_str} ---\n")
        except IOError as e:
            print(f"Warning: Could not open log file {log_path}: {e}. Logging disabled.", file=sys.stderr)
            stdout_handle = subprocess.PIPE
            stderr_handle = subprocess.PIPE
            log_file = None
    else:
        stdout_handle = subprocess.PIPE
        stderr_handle = subprocess.PIPE

    try:
        # Use shell=True ONLY if the command string requires it (like redirection)
        process = subprocess.run(cmd_str if shell else cmd_list,
                                 check=True, text=True,
                                 stdout=stdout_handle, stderr=stderr_handle,
                                 cwd=cwd, shell=shell) # Pass shell argument
        if log_file:
            stdout_handle.close()
        return process.returncode
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}: {cmd_str}", file=sys.stderr)
        if not log_file:
            if e.stderr: print(f"Stderr:\n{e.stderr}", file=sys.stderr)
            if e.stdout: print(f"Stdout:\n{e.stdout}", file=sys.stderr)
        if log_file: stdout_handle.close()
        return e.returncode
    except FileNotFoundError:
        # Extract the actual command if using shell
        cmd_to_report = cmd_str.split()[0] if shell else cmd_list[0]
        print(f"Error: Command not found: {cmd_to_report}", file=sys.stderr)
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f"\n--- ERROR: Command not found: {cmd_to_report} ---\n")
            stdout_handle.close()
        return 1
    except Exception as e:
        print(f"Error running command {cmd_str}: {e}", file=sys.stderr)
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f"\n--- ERROR: Exception running command: {e} ---\n")
            stdout_handle.close()
        return 1

def parse_snp_dists_matrix(matrix_file: Path) -> tuple[dict, list]:
    """
    Parses the snp-dists matrix output file.
    Assumes a tab-separated format where the first row is headers
    and the first column is also headers (matching the row).
    Returns a nested dictionary (hash_all) and a list of headers (head_all).
    Distances are expected to be proportions (0-1).
    """
    distances_hash = {}
    headers = []
    print(f"Parsing snp-dists matrix: {matrix_file}", file=sys.stderr)
    try:
        # Use pandas for robust parsing
        df = pd.read_csv(matrix_file, sep='\t', index_col=0) # Use first column as index
        headers = df.columns.tolist()
        if list(df.index) != headers:
             print(f"Warning: Row headers and Column headers in {matrix_file.name} do not match perfectly. Using column headers.", file=sys.stderr)
             # Potentially add more robust checking here if needed

        # Convert DataFrame to nested dictionary
        distances_hash = df.astype(float).to_dict('index')
        # Ensure symmetrical entries (snp-dists output should already be)
        for h1 in headers:
            for h2 in headers:
                # Ensure lookup exists both ways, handle potential NaN/missing if format is odd
                dist = distances_hash.get(h1, {}).get(h2)
                if dist is not None and np.isfinite(dist):
                     if h2 not in distances_hash: distances_hash[h2] = {}
                     distances_hash[h2][h1] = dist
                elif h1 != h2: # Don't warn for diagonal
                    print(f"Warning: Missing or non-finite distance between {h1} and {h2} in snp-dists matrix.", file=sys.stderr)
                    # Decide how to handle missing - set to max distance (1.0)? Skip?
                    # Setting to 1.0 might be safer for downstream calcs than NaN/Inf
                    if h1 in distances_hash: distances_hash[h1][h2] = 1.0
                    if h2 in distances_hash: distances_hash[h2][h1] = 1.0


    except pd.errors.EmptyDataError:
        print(f"Error: snp-dists matrix file {matrix_file.name} is empty.", file=sys.stderr)
        return {}, []
    except FileNotFoundError:
        print(f"Error: snp-dists matrix file {matrix_file.name} not found.", file=sys.stderr)
        return {}, []
    except Exception as e:
        print(f"Error parsing snp-dists matrix file {matrix_file.name}: {e}", file=sys.stderr)
        raise # Re-raise unexpected errors

    # Basic validation
    if not headers:
        print(f"Warning: No headers found in snp-dists matrix {matrix_file.name}.", file=sys.stderr)
    if not distances_hash:
        print(f"Warning: No distances parsed from snp-dists matrix {matrix_file.name}.", file=sys.stderr)

    return distances_hash, headers


def main():
    parser = argparse.ArgumentParser(description="HCV Transmission Test using snp-dists, Z-scores, and t-SNE plot.")
    parser.add_argument("fasta1_gz", type=Path, help="Path to the first gzipped FASTA file (sample 1).")
    parser.add_argument("fasta2_gz", type=Path, help="Path to the second gzipped FASTA file (sample 2).")
    parser.add_argument("base_dir", type=Path, help="Base directory for the pipeline (needed for Reports path).")
    parser.add_argument("tmp_dir", type=Path, help="Path to the shared temporary directory for this run.")
    parser.add_argument("--snp_dists_path", type=str, default="snp-dists", help="Path to the snp-dists executable.")
    parser.add_argument("--rscript_path", type=str, default="Rscript", help="Path to the Rscript executable.")
    parser.add_argument("--z_cutoff", type=float, default=0.8416, help="Z-score cutoff for classification (80th percentile).")
    parser.add_argument("--keep_tmp", action='store_true', help="Keep temporary run directory.")

    args = parser.parse_args()

    # --- Validate Inputs ---
    if not args.fasta1_gz.is_file():
        print(f"Error: Input FASTA 1 not found: {args.fasta1_gz}", file=sys.stderr); sys.exit(1)
    if not args.fasta2_gz.is_file():
        print(f"Error: Input FASTA 2 not found: {args.fasta2_gz}", file=sys.stderr); sys.exit(1)
    if not args.base_dir.is_dir():
        print(f"Error: Base directory not found: {args.base_dir}", file=sys.stderr); sys.exit(1)
    if not args.tmp_dir.is_dir():
        print(f"Error: Temporary directory not found: {args.tmp_dir}", file=sys.stderr); sys.exit(1)

    reports_dir = args.base_dir / "Reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    s1_name = args.fasta1_gz.name.split('.')[0]
    s2_name = args.fasta2_gz.name.split('.')[0]
    headers_s1 = set()
    headers_s2 = set()
    try:
        with gzip.open(args.fasta1_gz, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                headers_s1.add(record.id.replace('_R_', ''))
        with gzip.open(args.fasta2_gz, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                headers_s2.add(record.id.replace('_R_', ''))
    except Exception as e:
        print(f"Error reading headers from input FASTA files: {e}", file=sys.stderr); sys.exit(1)

    if not headers_s1: print(f"Warning: No sequences found in {args.fasta1_gz.name}", file=sys.stderr)
    if not headers_s2: print(f"Warning: No sequences found in {args.fasta2_gz.name}", file=sys.stderr)

    # --- Use Provided Temporary Directory ---
    temp_dir = args.tmp_dir
    print(f"Using temporary directory: {temp_dir}", file=sys.stderr)
    log_file = temp_dir / "run_log.txt"

    # --- Preprocessing: Unzip and Align ---
    print("Preprocessing: Unzipping and Aligning...", file=sys.stderr)
    tmp1_fa = temp_dir / "tmp1.fa"
    tmp2_fa = temp_dir / "tmp2.fa"
    data3_fasta = temp_dir / "data3.fasta"  # Combined unaligned
    data3_out = temp_dir / "data3.out"      # Combined aligned

    try:
        with gzip.open(args.fasta1_gz, "rb") as f_in, open(tmp1_fa, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        with gzip.open(args.fasta2_gz, "rb") as f_in, open(tmp2_fa, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        print(f"Error unzipping input files: {e}", file=sys.stderr); sys.exit(1)

    with open(data3_fasta, 'wb') as f_out, open(tmp1_fa, 'rb') as f_in1, open(tmp2_fa, 'rb') as f_in2:
        shutil.copyfileobj(f_in1, f_out)
        shutil.copyfileobj(f_in2, f_out)

    # Run MAFFT using Docker
    mafft_docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{temp_dir}:/data",
        "pegi3s/mafft",
        "mafft",
        "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder",
        f"/data/{data3_fasta.name}"
    ]

    print(f"Running MAFFT in Docker...", file=sys.stderr)
    try:
        with open(data3_out, 'w') as f_out_aln3, open(log_file, 'a') as flog:
            flog.write(f"\n--- Running MAFFT in Docker on {data3_fasta.name} ---\n")
            proc3 = subprocess.run(mafft_docker_cmd, cwd=temp_dir, text=True, capture_output=True, check=True)
            f_out_aln3.write(proc3.stdout)
            if proc3.stderr: flog.write(f"MAFFT Stderr:\n{proc3.stderr}\n")
        mafft_ret = 0
    except subprocess.CalledProcessError as e:
        print(f"MAFFT in Docker failed with exit code {e.returncode}. Check log file: {log_file}", file=sys.stderr)
        with open(log_file, "a") as flog: flog.write(f"MAFFT Error Output:\n{e.stderr}\n")
        mafft_ret = e.returncode
    except FileNotFoundError:
        print(f"Error: Docker not found.", file=sys.stderr); mafft_ret = 1
    except Exception as e:
        print(f"An unexpected error occurred running MAFFT in Docker: {e}", file=sys.stderr); mafft_ret = 1

    if mafft_ret != 0 or not data3_out.is_file() or data3_out.stat().st_size == 0:
        print("Error: MAFFT failed or produced an empty alignment file. Cannot calculate distances.", file=sys.stderr)
        print(f"\n{s1_name} vs {s2_name}:\tFailed to confirm data transmission (MAFFT alignment failed)")
        sys.exit(1)

    # --- Make Headers Unique Before Distance Calculation ---
    print("Making headers unique in combined alignment file...", file=sys.stderr)
    data3_unique_headers_out = temp_dir / "data3_unique_headers.out"
    header_count = 0
    try:
        with open(data3_out, "r") as f_in, open(data3_unique_headers_out, "w") as f_out:
            for record in SeqIO.parse(f_in, "fasta"):
                original_id_cleaned = record.id.replace('_R_', '')
                prefix = ""
                if original_id_cleaned in headers_s1:
                    prefix = s1_name
                elif original_id_cleaned in headers_s2:
                    prefix = s2_name
                else:
                    print(f"Warning: Header {record.id} not found in original sample lists.", file=sys.stderr)
                    prefix = "UNKNOWN"

                record.id = f"{prefix}_{record.id}" # Prepend sample name
                record.description = "" # Clear description
                SeqIO.write(record, f_out, "fasta")
                header_count += 1
        if header_count < 2:
             print(f"Warning: Combined alignment {data3_unique_headers_out.name} has fewer than 2 sequences after header processing.", file=sys.stderr)
             # Potentially exit early if 0 or 1 sequences remain

    except Exception as e:
        print(f"Error creating unique headers file: {e}", file=sys.stderr); sys.exit(1)

    # --- Calculate Pairwise Distances using snp-dists ---
    print("Calculating pairwise distances using snp-dists...", file=sys.stderr)
    pmatrix_file = temp_dir / "ident.pmatrix" # Output file for snp-dists
    mmatrix_file = temp_dir / "ident.mmatrix" # Output file for snp-dists
    # Check if snp-dists path needs quoting (e.g., if it contains spaces)
    # Assuming simple path or in PATH for now. Add quoting if needed.
    # Use shell=True because of the '>' redirection
    snp_dists_cmd = f"{args.snp_dists_path} -b {data3_unique_headers_out.resolve()} > {pmatrix_file.resolve()}"
    snp_ret = run_command(snp_dists_cmd, cwd=temp_dir, log_file=log_file, shell=True)
    snp_dists_cmd2 = f"{args.snp_dists_path} -m {data3_unique_headers_out.resolve()} > {mmatrix_file.resolve()}"
    run_command(snp_dists_cmd2, cwd=temp_dir, log_file=log_file, shell=True)
    hash_all = {}
    head_all = []
    if snp_ret != 0:
        print(f"Error: snp-dists failed (exit code {snp_ret}). Cannot proceed.", file=sys.stderr)
        print(f"\n{s1_name} vs {s2_name}:\tFailed to confirm data transmission (snp-dists failed)")
        sys.exit(1)
    elif not pmatrix_file.is_file() or pmatrix_file.stat().st_size == 0:
        print(f"Error: snp-dists ran but output matrix {pmatrix_file.name} is missing or empty.", file=sys.stderr)
        print(f"\n{s1_name} vs {s2_name}:\tFailed to confirm data transmission (snp-dists matrix empty)")
        sys.exit(1)
    else:
        # --- Parse the snp-dists output matrix ---
        try:
            hash_all, head_all = parse_snp_dists_matrix(pmatrix_file)
        except Exception as e:
             print(f"Fatal Error: Failed to parse snp-dists output. {e}", file=sys.stderr)
             print(f"\n{s1_name} vs {s2_name}:\tFailed to confirm data transmission (matrix parsing failed)")
             sys.exit(1)

    # Validate parsed results
    if not head_all or not hash_all:
         print("Error: Failed to parse headers or distances from snp-dists output.", file=sys.stderr)
         print(f"\n{s1_name} vs {s2_name}:\tFailed to confirm data transmission (matrix parsing failed)")
         sys.exit(1)
    if len(head_all) < 2:
         print(f"Warning: Fewer than 2 sequences found in snp-dists matrix ({len(head_all)}). Z-score calculations might be unreliable.", file=sys.stderr)
         # Don't exit here, let downstream handle few sequences, but stats might be zero


    # --- Filter distances for within-sample calculations ---
    print("Filtering within-sample distances...", file=sys.stderr)
    array_d1 = []
    array_d2 = []

    # Use the headers and hash directly from snp-dists output
    num_seqs_from_snp_dists = len(head_all)
    for i in range(num_seqs_from_snp_dists):
        h1 = head_all[i]
        for j in range(i + 1, num_seqs_from_snp_dists):
            h2 = head_all[j]
            # Distance lookup from the parsed hash
            dist = hash_all.get(h1, {}).get(h2)

            if dist is None:
                print(f"Warning: Missing distance between {h1} and {h2} in parsed snp-dists hash", file=sys.stderr)
                dist = 1.0 # Assign max distance if missing? Or skip?

            # Check prefixes for within-sample calculation
            h1_prefix = h1.split('_', 1)[0]
            h2_prefix = h2.split('_', 1)[0]
            if h1_prefix == s1_name and h2_prefix == s1_name:
                array_d1.append(dist)
            elif h1_prefix == s2_name and h2_prefix == s2_name:
                array_d2.append(dist)

    # Check if any within-sample distances were found
    # Compare unique sequence count in original sample vs effective count in alignment
    # Get effective counts after alignment/snp-dists
    unique_s1_in_final = sum(1 for h in head_all if h.startswith(s1_name + '_'))
    unique_s2_in_final = sum(1 for h in head_all if h.startswith(s2_name + '_'))

    if not array_d1 and unique_s1_in_final > 1:
        print(f"Warning: No within-sample distances calculated for {s1_name} (unique headers in matrix: {unique_s1_in_final}). Check alignment/snp-dists output.", file=sys.stderr)
    if not array_d2 and unique_s2_in_final > 1:
        print(f"Warning: No within-sample distances calculated for {s2_name} (unique headers in matrix: {unique_s2_in_final}). Check alignment/snp-dists output.", file=sys.stderr)


    # --- Calculate Within-Sample Stats ---
    # Ensure array is numpy array for calculations
    np_array_d1 = np.array(array_d1, dtype=float)
    np_array_d2 = np.array(array_d2, dtype=float)

    mean_d1 = np.mean(np_array_d1) if np_array_d1.size > 0 else 0.0
    stdev_d1 = np.std(np_array_d1) if np_array_d1.size > 0 else 0.0
    mean_d2 = np.mean(np_array_d2) if np_array_d2.size > 0 else 0.0
    stdev_d2 = np.std(np_array_d2) if np_array_d2.size > 0 else 0.0

    # Prevent division by zero if standard deviation is very small or zero
    stdev_d1 = max(stdev_d1, 1e-9) # Use max instead of conditional assignment
    stdev_d2 = max(stdev_d2, 1e-9)


    # --- Classify Sequences & Write Stats File ---
    stats_file_path = temp_dir / "stats.txt"
    print(f"Classifying sequences based on Z-scores and writing stats to: {stats_file_path}", file=sys.stderr)
    results = {s1_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0},
               s2_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0}}
    legend_data = [] # Still needed for R plot generation later IF using internal legend gen

    try:
        with open(stats_file_path, "w") as f_stats:
            f_stats.write("--- Intra-Sample Statistics ---\n")
            f_stats.write(f"Sample\tMean_Distance (from snp-dists)\tStdev_Distance\n")
            f_stats.write(f"{s1_name}\t{'N/A' if mean_d1 == 0 else f'{mean_d1:.6f}'}\t{'N/A' if stdev_d1 <= 1e-9 else f'{stdev_d1:.6f}'}\n")
            f_stats.write(f"{s2_name}\t{'N/A' if mean_d2 == 0 else f'{mean_d2:.6f}'}\t{'N/A' if stdev_d2 <= 1e-9 else f'{stdev_d2:.6f}'}\n")
            f_stats.write("\n--- Per-Sequence Z-Scores ---\n")
            f_stats.write("Sequence_Header\tMean_Z_vs_Sample1\tMean_Z_vs_Sample2\n")

            # Use head_all from snp-dists parsing
            for h_current in head_all:
                current_sample_name = None
                original_header_part = None
                legend_cat = 0 # Placeholder

                if h_current.startswith(s1_name + '_'):
                    current_sample_name = s1_name
                    original_header_part = h_current[len(s1_name)+1:]
                    legend_cat = 2 # Category for sample 1
                elif h_current.startswith(s2_name + '_'):
                    current_sample_name = s2_name
                    original_header_part = h_current[len(s2_name)+1:]
                    legend_cat = 4 # Category for sample 2
                else:
                    print(f"Warning: Could not determine original sample for header {h_current} from snp-dists matrix. Skipping classification.", file=sys.stderr)
                    continue # Skip sequences that don't match expected prefixes

                # --- Extract count for size calculation (IF NEEDED for plot legend) ---
                # This assumes count is reliably the last part after '_'
                try:
                    count_str = h_current.split('_')[-1]
                    # Further check if it looks like a typical sequence ID vs. just a count
                    # Let's assume simple count extraction for now
                    count = int(count_str)
                    scaled_size = math.log10(count) if count > 0 else 0
                except (IndexError, ValueError):
                    print(f"Warning: Could not extract numeric count from unique header {h_current}. Using default size 1 for plot legend.", file=sys.stderr)
                    scaled_size = 1 # Default size

                legend_data.append((h_current, legend_cat, scaled_size))
                # ---- End Count Extraction ----

                z_scores_vs_s1 = []
                z_scores_vs_s2 = []

                # Compare h_current to all OTHER sequences using snp-dists distances
                for h_other in head_all:
                    if h_current == h_other: continue # Skip self-comparison

                    dist = hash_all.get(h_current, {}).get(h_other)
                    if dist is None:
                       print(f"Warning: Missing distance between {h_current} and {h_other} during Z-score calculation (should not happen if parsing was complete).", file=sys.stderr)
                       continue # Skip if somehow missing

                    # Accumulate scores based on the sample h_other belongs to
                    if h_other.startswith(s1_name + '_'):
                        # Only calculate Z if mean/stdev for Sample 1 are valid
                        if unique_s1_in_final > 1:
                             z_scores_vs_s1.append((dist - mean_d1) / stdev_d1)
                    elif h_other.startswith(s2_name + '_'):
                        # Only calculate Z if mean/stdev for Sample 2 are valid
                        if unique_s2_in_final > 1:
                             z_scores_vs_s2.append((dist - mean_d2) / stdev_d2)

                # Calculate mean Z-scores
                mean_z_vs_s1 = np.mean(z_scores_vs_s1) if z_scores_vs_s1 else float('inf')
                mean_z_vs_s2 = np.mean(z_scores_vs_s2) if z_scores_vs_s2 else float('inf')

                # Handle cases where no comparisons could be made (e.g., only one sequence from a sample)
                if not z_scores_vs_s1: mean_z_vs_s1 = float('inf')
                if not z_scores_vs_s2: mean_z_vs_s2 = float('inf')


                z1_str = f"{mean_z_vs_s1:.6f}" if np.isfinite(mean_z_vs_s1) else "Inf"
                z2_str = f"{mean_z_vs_s2:.6f}" if np.isfinite(mean_z_vs_s2) else "Inf"
                f_stats.write(f"{h_current}\t{z1_str}\t{z2_str}\n")


                # --- Classification Logic (same as before, uses mean_z_vs_s1/s2) ---
                if current_sample_name == s1_name:
                    q1 = 1 if not np.isfinite(mean_z_vs_s1) or mean_z_vs_s1 >= args.z_cutoff else 0
                    q2 = 1 if not np.isfinite(mean_z_vs_s2) or mean_z_vs_s2 >= args.z_cutoff else 0
                    if q1 == 0 and q2 == 0: results[s1_name]['AA'] += 1
                    elif q1 == 1 and q2 == 0: results[s1_name]['BA'] += 1
                    elif q1 == 0 and q2 == 1: results[s1_name]['AB'] += 1
                    else: results[s1_name]['BB'] += 1 # q1=1, q2=1 (or Inf)
                elif current_sample_name == s2_name:
                    q2_perl = 1 if not np.isfinite(mean_z_vs_s1) or mean_z_vs_s1 >= args.z_cutoff else 0 # vs S1
                    q1_perl = 1 if not np.isfinite(mean_z_vs_s2) or mean_z_vs_s2 >= args.z_cutoff else 0 # vs S2
                    if q1_perl == 0 and q2_perl == 0: results[s2_name]['AA'] += 1
                    elif q1_perl == 1 and q2_perl == 0: results[s2_name]['AB'] += 1
                    elif q1_perl == 0 and q2_perl == 1: results[s2_name]['BA'] += 1
                    else: results[s2_name]['BB'] += 1 # q1=1, q2=1 (or Inf)
                # --- End Classification Logic ---

    except IOError as e:
        print(f"Error writing stats file {stats_file_path}: {e}", file=sys.stderr)
    except Exception as e:
        print(f"Error during classification or stats file writing: {e}", file=sys.stderr)
        # Ensure results dict exists even if loop failed partially
        if 'results' not in locals():
             results = {s1_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0},
                        s2_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0}}

    # --- Determine Final Results ---
    print("Determining final classification...", file=sys.stderr)
    # (Classification logic remains the same as previous version)
    res1 = results[s1_name]
    tot1 = sum(res1.values())
    result_s1 = "Failed"
    # Add check for tot1 > 0 before division
    if tot1 > 0:
        # Adjust thresholds if necessary based on snp-dists results behavior
        if res1['AA'] / tot1 >= 0.8: result_s1 = "Mixed"
        elif res1['AB'] / tot1 >= 0.5: result_s1 = "Homogenious" # Note: Typo in original script 'Homogenious' kept for consistency
        elif res1['BA'] / tot1 >= 0.5: result_s1 = "recipient"
        elif res1['BB'] / tot1 > 0.5 : # Check if BB is the majority if others fail
              result_s1 = "Too many unclassified reads"
        # Add fallback if no category meets criteria?
        elif result_s1 == "Failed": # If still failed after checks
             print(f"Warning: Could not classify sample {s1_name} based on thresholds.", file=sys.stderr)
             # Optionally assign a default or leave as Failed


    res2 = results[s2_name]
    tot2 = sum(res2.values())
    result_s2 = "Failed"
    if tot2 > 0:
        # Adjust thresholds if necessary
        if res2['AA'] / tot2 >= 0.2: result_s2 = "Mixed" # Note the lower threshold here compared to S1
        elif res2['AB'] / tot2 >= 0.5: result_s2 = "recipient"
        elif res2['BA'] / tot2 >= 0.5: result_s2 = "Homogenious" # Typo kept
        elif res2['BB'] / tot2 > 0.5 : # Check if BB is the majority
              result_s2 = "Too many unclassified reads"
        elif result_s2 == "Failed":
             print(f"Warning: Could not classify sample {s2_name} based on thresholds.", file=sys.stderr)

    # --- Print Summary to STDOUT ---
    print("\n--- Classification Summary ---")
    print("Sample\tIntermixture(AA)\tHomog(AB/BA)\tRecip(BA/AB)\tUnclustered(BB)\tTotal")
    print("----------------------------------------------------------------------------")
    # Clarify meaning based on original sample origin
    print(f"{s1_name}:\t{res1['AA']}(Mixed)\t\t{res1['AB']}(Homog S1)\t{res1['BA']}(Like S2)\t{res1['BB']}\t\t{tot1}")
    print(f"{s2_name}:\t{res2['AA']}(Mixed)\t\t{res2['AB']}(Like S1)\t{res2['BA']}(Homog S2)\t{res2['BB']}\t\t{tot2}")
    print(f"\nIntermediate Classifications: {result_s1} / {result_s2}\n")

    # --- Determine Final Conclusion ---
    # (Logic remains the same)
    conclusion = "Could Not Determine"
    # Check for explicit failure conditions first
    if "Failed" in result_s1 or "Failed" in result_s2 or "unclassified" in result_s1 or "unclassified" in result_s2:
       conclusion = "Failed to confirm data transmission (classification issues)"
    elif result_s1 == "Mixed" and result_s2 == "Mixed": conclusion = "Complete intersample mixture"
    elif result_s1 == "Mixed" and result_s2 == "Homogenious": conclusion = f"Transmission likely from {s2_name} to {s1_name}"
    elif result_s2 == "Mixed" and result_s1 == "Homogenious": conclusion = f"Transmission likely from {s1_name} to {s2_name}"
    elif result_s1 == "Homogenious" and result_s2 == "Homogenious": conclusion = "No transmission detected (distinct populations)"
    # Interpret 'recipient' in context
    elif result_s1 == "recipient" and result_s2 == "Homogenious": conclusion = f"Transmission likely from {s2_name} to {s1_name}" # S1 looks like S2, S2 is distinct
    elif result_s2 == "recipient" and result_s1 == "Homogenious": conclusion = f"Transmission likely from {s1_name} to {s2_name}" # S2 looks like S1, S1 is distinct
    # Add cases for other combinations if needed, e.g., Mixed/Recipient
    elif result_s1 == "Mixed" and result_s2 == "recipient": conclusion = f"Transmission likely from {s1_name} to {s2_name} (S2 has S1 background)"
    elif result_s2 == "Mixed" and result_s1 == "recipient": conclusion = f"Transmission likely from {s2_name} to {s1_name} (S1 has S2 background)"
    elif result_s1 == "recipient" and result_s2 == "recipient": conclusion = "Complex relationship or bidirectional transmission possible"


    print(f"--- Final Conclusion ---")
    print(f"{s1_name} vs {s2_name}:")
    print(conclusion)
    print("______________________________________\n")

    # --- Generate t-SNE plot using External R Script ---
    # The snp-dists matrix (pmatrix_file) was already generated earlier.
    # We just need to run the R script now.
    print("Generating t-SNE plot using external R script...", file=sys.stderr)
    tsne_input_matrix_file = mmatrix_file # Use the already generated matrix
    # Make external_r_script_path configurable or find relative to script? Hardcoded for now.
    external_r_script_path = Path("/Volumes/MacOS_Storage/Projects/HCV_pipeline/plot_tsne.R")
    # Check if R script exists
    if not external_r_script_path.is_file():
        print(f"Error: R script not found at {external_r_script_path}. Cannot generate plot.", file=sys.stderr)
    else:
        temp_pdf_output_name = "ident_tsne.pdf" # Expected output name from plot_tsne.R
        temp_pdf_output_path = temp_dir / temp_pdf_output_name
        final_pdf_output_path = reports_dir / f"Rtsne_{s1_name}_{s2_name}.pdf"

        print(f"Executing R script: {external_r_script_path} with input {tsne_input_matrix_file.name}", file=sys.stderr)
        # Use Rscript path from arguments
        r_ret = run_command([args.rscript_path, str(external_r_script_path), str(tsne_input_matrix_file)],
                            cwd=temp_dir, log_file=log_file) # Run R script without shell


        # Check R script result and move/rename the PDF
        if r_ret == 0 and temp_pdf_output_path.is_file():
            print(f"R script finished. Moving plot from {temp_pdf_output_path} to {final_pdf_output_path}", file=sys.stderr)
            try:
                shutil.move(str(temp_pdf_output_path), str(final_pdf_output_path))
                print(f"Plot successfully saved to: {final_pdf_output_path}", file=sys.stderr)
            except Exception as e:
                print(f"Error moving plot {temp_pdf_output_path} to {final_pdf_output_path}: {e}", file=sys.stderr)
        elif r_ret != 0:
            print(f"Warning: External R script execution failed (exit code {r_ret}). Plot may not be generated. Check log file: {log_file}", file=sys.stderr)
        else: # r_ret == 0 but file doesn't exist
            print(f"Warning: R script finished successfully, but expected output PDF not found: {temp_pdf_output_path}", file=sys.stderr)

    # Cleanup is handled by the calling script via args.keep_tmp

if __name__ == "__main__":
    # Add pandas as a dependency requirement
    try:
        import pandas
    except ImportError:
        print("Error: This script requires the 'pandas' library. Please install it (e.g., 'pip install pandas')", file=sys.stderr)
        sys.exit(1)
    main()