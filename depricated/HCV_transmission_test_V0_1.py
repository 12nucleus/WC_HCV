# HCV_transmission_test_V0_1.py
# Optimized version using NumPy for distance calculation
import argparse
import subprocess
import gzip
import shutil
import tempfile
from pathlib import Path
import sys
import os
import math
import numpy as np  # Import NumPy
from Bio import SeqIO, AlignIO


def calculate_pairwise_distances_optimized(aligned_fasta_path: Path) -> tuple[dict, list, list]:
    """
    Reads an aligned FASTA, calculates pairwise distances using the custom metric
    (optimized with NumPy), and returns the distance hash, list of distances,
    and list of headers.
    """
    try:
        # Read alignment - consider memory for very large alignments
        alignment = AlignIO.read(aligned_fasta_path, "fasta")
        num_seqs = len(alignment)
        if num_seqs == 0:
            print(f"Warning: Alignment {aligned_fasta_path.name} is empty.", file=sys.stderr)
            return {}, [], []
        aln_len = alignment.get_alignment_length()
    except ValueError as e:
        print(f"Warning: Could not read alignment {aligned_fasta_path}. Might be empty or invalid format. {e}", file=sys.stderr)
        if "Empty file" in str(e) or "No records found" in str(e):
            return {}, [], []
        else:
            raise
    except Exception as e:
        print(f"Error reading alignment file {aligned_fasta_path}: {e}", file=sys.stderr)
        raise

    headers = [record.id.replace('_R_', '') for record in alignment]  # Clean header like Perl script

    # Convert alignment to NumPy array of bytes for efficient comparison
    # This step can consume memory proportional to N * L
    try:
        aln_array = np.array([list(str(record.seq).upper().encode('ascii')) for record in alignment], dtype='S1')
    except UnicodeEncodeError:
        print(f"Error: Non-ASCII characters found in alignment {aligned_fasta_path.name}. Cannot process with NumPy 'S1' dtype.", file=sys.stderr)
        # Fallback or different dtype might be needed if non-ASCII is expected
        raise

    distances_list = []
    distances_hash = {h: {} for h in headers}

    if num_seqs < 2:
        print(f"Warning: Alignment {aligned_fasta_path.name} has less than 2 sequences. No pairwise distances calculated.", file=sys.stderr)
        return distances_hash, distances_list, headers

    # Pre-calculate N positions and valid character positions for efficiency
    is_n = (aln_array == b'N')
    valid_chars = np.array([b'A', b'T', b'C', b'G', b'-'], dtype='S1')
    is_valid = np.isin(aln_array, valid_chars)

    # Iterate through unique pairs (i, j) where j > i
    for i in range(num_seqs):
        seq_i = aln_array[i, :]
        is_n_i = is_n[i, :]
        is_valid_i = is_valid[i, :]
        n1 = np.sum(is_n_i)  # Count Ns in sequence i

        for j in range(i + 1, num_seqs):
            seq_j = aln_array[j, :]
            is_n_j = is_n[j, :]
            is_valid_j = is_valid[j, :]
            h1 = headers[i]
            h2 = headers[j]

            # Calculate terms using vectorized NumPy operations
            mismatches = (seq_i != seq_j)
            both_n = is_n_i & is_n_j
            both_valid = is_valid_i & is_valid_j

            same = np.sum(mismatches)  # hd equivalent
            nn = np.sum(both_n)  # hdn equivalent
            # n1 is precalculated
            n2 = np.sum(is_n_j)  # n equivalent for seq j
            v = np.sum(both_valid)  # valid equivalent

            # Calculate 'end' term from Perl script
            end = same - ((n1 - nn) + (n2 - nn))

            # Calculate final distance 'pp'
            denominator = end + v
            if denominator == 0:
                pp = 0.0
            elif denominator < 0:
                # This shouldn't happen with the formula, but add safety check
                print(f"Warning: Negative denominator ({denominator}) calculating distance between {h1} and {h2}. Setting distance to 1.0.", file=sys.stderr)
                pp = 1.0
            else:
                pp = end / denominator

            # Clamp distance between 0 and 1
            pp = max(0.0, min(1.0, pp))

            distances_list.append(pp)
            distances_hash[h1][h2] = pp
            distances_hash[h2][h1] = pp

    return distances_hash, distances_list, headers


def run_command(cmd_list, cwd=None, log_file=None):
    """Runs an external command, logs output, checks return code."""
    print(f"Running command: {' '.join(cmd_list)}", file=sys.stderr)
    stdout_pipe = subprocess.PIPE
    stderr_pipe = subprocess.PIPE
    if log_file:
        # Ensure log_file is Path object or string
        log_path = Path(log_file)
        try:
            stdout_handle = open(log_path, 'a')  # Append stdout
            stderr_handle = stdout_handle  # Append stderr to the same file
            stdout_handle.write(f"\n--- Running: {' '.join(cmd_list)} ---\n")  # Add separator
        except IOError as e:
            print(f"Warning: Could not open log file {log_path}: {e}. Logging disabled for this command.", file=sys.stderr)
            stdout_handle = subprocess.PIPE
            stderr_handle = subprocess.PIPE
            log_file = None  # Disable closing handle later
    else:
        stdout_handle = subprocess.PIPE
        stderr_handle = subprocess.PIPE

    try:
        process = subprocess.run(cmd_list, check=True, text=True,
                                 stdout=stdout_handle, stderr=stderr_handle, cwd=cwd)
        if log_file:
            stdout_handle.close()
        return process.returncode
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}: {' '.join(cmd_list)}", file=sys.stderr)
        # If logging to file, error details are already appended by run()
        if not log_file:  # Only print if not logging to file
            if e.stderr:
                print(f"Stderr:\n{e.stderr}", file=sys.stderr)
            if e.stdout:
                print(f"Stdout:\n{e.stdout}", file=sys.stderr)
        if log_file:
            stdout_handle.close()
        return e.returncode
    except FileNotFoundError:
        print(f"Error: Command not found: {cmd_list[0]}", file=sys.stderr)
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f"\n--- ERROR: Command not found: {cmd_list[0]} ---\n")
            stdout_handle.close()
        return 1
    except Exception as e:
        print(f"Error running command {' '.join(cmd_list)}: {e}", file=sys.stderr)
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f"\n--- ERROR: Exception running command: {e} ---\n")
            stdout_handle.close()
        return 1


def main():
    parser = argparse.ArgumentParser(description="HCV Transmission Test using sequence alignment, distance, and Z-scores (Optimized).")
    parser.add_argument("fasta1_gz", type=Path, help="Path to the first gzipped FASTA file (sample 1).")
    parser.add_argument("fasta2_gz", type=Path, help="Path to the second gzipped FASTA file (sample 2).")
    parser.add_argument("base_dir", type=Path, help="Base directory for the pipeline (needed for Reports path).")
    parser.add_argument("tmp_dir", type=Path, help="Path to the shared temporary directory for this run.") # New argument
    parser.add_argument("--snp_dists_path", type=str, default="snp-dists", help="Path to the snp-dists executable.")
    parser.add_argument("--rscript_path", type=str, default="Rscript", help="Path to the Rscript executable.")
    parser.add_argument("--z_cutoff", type=float, default=0.8416, help="Z-score cutoff for classification (80th percentile).")
    parser.add_argument("--keep_tmp", action='store_true', help="Keep temporary run directory.")

    args = parser.parse_args()

    # --- Validate Inputs ---
    if not args.fasta1_gz.is_file():
        print(f"Error: Input FASTA 1 not found: {args.fasta1_gz}", file=sys.stderr)
        sys.exit(1)
    if not args.fasta2_gz.is_file():
        print(f"Error: Input FASTA 2 not found: {args.fasta2_gz}", file=sys.stderr)
        sys.exit(1)
    if not args.base_dir.is_dir(): # Keep check for base_dir needed for Reports
        print(f"Error: Base directory not found: {args.base_dir}", file=sys.stderr)
        sys.exit(1)
    if not args.tmp_dir.is_dir(): # Check the provided tmp_dir exists
        print(f"Error: Temporary directory not found: {args.tmp_dir}", file=sys.stderr)
        sys.exit(1)

    reports_dir = args.base_dir / "Reports" # Reports still go relative to base_dir
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Extract sample names and read headers *before* temp dir creation
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
        print(f"Error reading headers from input FASTA files: {e}", file=sys.stderr)
        sys.exit(1)

    if not headers_s1:
        print(f"Warning: No sequences found in {args.fasta1_gz.name}", file=sys.stderr)
        # Decide if this is fatal or should continue with empty set
    if not headers_s2:
        print(f"Warning: No sequences found in {args.fasta2_gz.name}", file=sys.stderr)

    # --- Use Provided Temporary Directory ---
    # No longer create a temp dir here, use the one provided.
    # try: # Remove the try block associated with temp_dir_manager
    temp_dir = args.tmp_dir # Use the passed-in directory
    print(f"Using temporary directory: {temp_dir}", file=sys.stderr)
    log_file = temp_dir / "run_log.txt"

    # --- Preprocessing: Unzip and Align ---
    print("Preprocessing: Unzipping and Aligning...", file=sys.stderr)
    tmp1_fa = temp_dir / "tmp1.fa" # Still needed for combining
    tmp2_fa = temp_dir / "tmp2.fa" # Still needed for combining
    # data1_out and data2_out definitions are removed as they are not created by MAFFT
    data3_fasta = temp_dir / "data3.fasta"  # Combined unaligned
    data3_out = temp_dir / "data3.out"  # Combined aligned (This is the only alignment output)

    try:
        with gzip.open(args.fasta1_gz, "rb") as f_in, open(tmp1_fa, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        with gzip.open(args.fasta2_gz, "rb") as f_in, open(tmp2_fa, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        print(f"Error unzipping input files: {e}", file=sys.stderr)
        sys.exit(1)

    with open(data3_fasta, 'wb') as f_out, open(tmp1_fa, 'rb') as f_in1, open(tmp2_fa, 'rb') as f_in2:
        shutil.copyfileobj(f_in1, f_out)
        shutil.copyfileobj(f_in2, f_out)

    # Run MAFFT using Docker
    mafft_docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{temp_dir}:/data",
        "pegi3s/mafft", # Image name
        "mafft",        # Explicitly specify the command to run inside the container
        "--adjustdirection", "--auto", "--quiet", "--thread", "1", "--reorder", # MAFFT options
        f"/data/{data3_fasta.name}" # Input file inside container
    ]

    print(f"Running MAFFT in Docker...", file=sys.stderr)
    try:
        with open(data3_out, 'w') as f_out_aln3, open(log_file, 'a') as flog:
            flog.write(f"\n--- Running MAFFT in Docker on {data3_fasta.name} ---\n")
            proc3 = subprocess.run(mafft_docker_cmd, cwd=temp_dir, text=True, capture_output=True, check=True)
            f_out_aln3.write(proc3.stdout)
            if proc3.stderr:
                flog.write(f"MAFFT Stderr:\n{proc3.stderr}\n")
        mafft_ret = 0
    except subprocess.CalledProcessError as e:
        print(f"MAFFT in Docker failed with exit code {e.returncode}. Check log file: {log_file}", file=sys.stderr)
        with open(log_file, "a") as flog:  # Log error details if possible
            flog.write(f"MAFFT Error Output:\n{e.stderr}\n")
        mafft_ret = e.returncode
    except FileNotFoundError:
        print(f"Error: Docker not found.", file=sys.stderr)
        mafft_ret = 1
    except Exception as e:
        print(f"An unexpected error occurred running MAFFT in Docker: {e}", file=sys.stderr)
        mafft_ret = 1

    # Check MAFFT result *before* trying to calculate distances
    if mafft_ret != 0 or not data3_out.is_file() or data3_out.stat().st_size == 0:
        print("Error: MAFFT failed or produced an empty alignment file. Cannot calculate distances.", file=sys.stderr)
        # Clean up and exit if MAFFT failed critically
        # temp_dir_manager.cleanup() # Master script handles cleanup
        print(f"\n{s1_name} vs {s2_name}:\tFailed to confirm data transmission (MAFFT alignment failed)")
        sys.exit(1) # Exit with error code

    # --- Make Headers Unique Before Distance Calculation ---
    print("Making headers unique in combined alignment file...", file=sys.stderr)
    data3_unique_headers_out = temp_dir / "data3_unique_headers.out"
    try:
        with open(data3_out, "r") as f_in, open(data3_unique_headers_out, "w") as f_out:
            for record in SeqIO.parse(f_in, "fasta"):
                original_id_cleaned = record.id.replace('_R_', '') # Match cleaning done earlier
                if original_id_cleaned in headers_s1:
                    record.id = f"{s1_name}_{record.id}"
                elif original_id_cleaned in headers_s2:
                    record.id = f"{s2_name}_{record.id}"
                else:
                    # Should not happen if headers_s1/s2 are correct, but add warning
                    print(f"Warning: Header {record.id} not found in original sample lists.", file=sys.stderr)
                    record.id = f"UNKNOWN_{record.id}"
                record.description = "" # Clear description to avoid extra info
                SeqIO.write(record, f_out, "fasta")
    except Exception as e:
        print(f"Error creating unique headers file: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Calculate Pairwise Distances (Optimized) ---
    print("Calculating pairwise distances from unique-header alignment...", file=sys.stderr)
    # Calculate all distances from the file with unique headers
    hash_all, array_all, head_all = calculate_pairwise_distances_optimized(data3_unique_headers_out)

    if not head_all:
         print("Error: Combined alignment resulted in no sequences after processing.", file=sys.stderr)
         # temp_dir_manager.cleanup() # Master script handles cleanup
         print(f"\n{s1_name} vs {s2_name}:\tFailed to confirm data transmission (empty combined alignment)")
         sys.exit(1) # Exit with error code

    # Filter distances for within-sample calculations
    print("Filtering within-sample distances...", file=sys.stderr)
    array_d1 = []
    array_d2 = []
    # head_map = {h: idx for idx, h in enumerate(head_all)} # For quick lookup if needed
    
    for i in range(len(head_all)):
        h1 = head_all[i]
        
        for j in range(i + 1, len(head_all)):
            h2 = head_all[j]
            # Distance should exist in hash_all if calculated
            dist = hash_all.get(h1, {}).get(h2)
            if dist is None:
                # This shouldn't happen if hash_all is complete, but good check
                print(f"Warning: Missing distance between {h1} and {h2} in hash_all", file=sys.stderr)
                continue

            # Check if both unique headers belong to the same original sample based on prefix
            h1_prefix = h1.split('_', 1)[0]
            h2_prefix = h2.split('_', 1)[0]
            if h1_prefix == s1_name and h2_prefix == s1_name:
                array_d1.append(dist)
            elif h1_prefix == s2_name and h2_prefix == s2_name:
                array_d2.append(dist)

    # Check if we got any within-sample distances (using the filtered arrays)
    if not array_d1 and len(headers_s1) > 1 :
         print(f"Warning: No within-sample distances calculated for {s1_name} (original headers: {len(headers_s1)}). Check alignment.", file=sys.stderr)
    if not array_d2 and len(headers_s2) > 1:
         print(f"Warning: No within-sample distances calculated for {s2_name} (original headers: {len(headers_s2)}). Check alignment.", file=sys.stderr)

    # --- Calculate Within-Sample Stats ---
    mean_d1 = np.mean(array_d1) if array_d1 else 0
    stdev_d1 = np.std(array_d1) if array_d1 else 0
    mean_d2 = np.mean(array_d2) if array_d2 else 0
    stdev_d2 = np.std(array_d2) if array_d2 else 0
    #print (f"xxxx array_d1: {array_d1}, mean_d1: {mean_d1}, stdev_d1: {stdev_d1}", file=sys.stderr)
    stdev_d1 = stdev_d1 if stdev_d1 > 1e-9 else 1.0
    stdev_d2 = stdev_d2 if stdev_d2 > 1e-9 else 1.0

    # --- Classify Sequences & Write Stats File ---
    stats_file_path = temp_dir / "stats.txt"
    print(f"Classifying sequences based on Z-scores and writing stats to: {stats_file_path}", file=sys.stderr)
    results = {s1_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0},
               s2_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0}}
    legend_data = []

    try:
        with open(stats_file_path, "w") as f_stats:
            # Write overall stats first
            f_stats.write("--- Intra-Sample Statistics ---\n")
            f_stats.write(f"Sample\tMean_Distance\tStdev_Distance\n")
            f_stats.write(f"{s1_name}\t{mean_d1:.6f}\t{stdev_d1:.6f}\n")
            f_stats.write(f"{s2_name}\t{mean_d2:.6f}\t{stdev_d2:.6f}\n")
            f_stats.write("\n--- Per-Sequence Z-Scores ---\n")
            f_stats.write("Sequence_Header\tMean_Z_vs_Sample1\tMean_Z_vs_Sample2\n")

            # Analyze sequences based on the unique headers from the final alignment (head_all)
            for h_current in head_all:
                # Determine the original sample for the current header
                current_sample_name = None
                original_header_part = None # Store the part after the sample prefix if needed
                if h_current.startswith(s1_name + '_'):
                    current_sample_name = s1_name
                    original_header_part = h_current[len(s1_name)+1:]
                    legend_cat = 2 # Category for sample 1
                elif h_current.startswith(s2_name + '_'):
                    current_sample_name = s2_name
                    original_header_part = h_current[len(s2_name)+1:]
                    legend_cat = 4 # Category for sample 2
                else:
                    # This case should ideally not happen if unique header creation worked
                    print(f"Warning: Could not determine original sample for unique header {h_current}. Skipping classification.", file=sys.stderr)
                    continue

                # Extract count for size calculation (using unique header)
                try:
                    # Extract last part after underscore, assuming it's the count
                    count = int(h_current.split('_')[-1])
                    # Use log10 for size as done previously for sample 1
                    scaled_size = math.log10(count) if count > 0 else 0
                except (IndexError, ValueError):
                    print(f"Warning: Could not extract count from unique header {h_current}. Using default size 1.", file=sys.stderr)
                    scaled_size = 1 # Default size if count extraction fails

                # Append data needed for potential plotting (using unique header)
                legend_data.append((h_current, legend_cat, scaled_size))

                # Calculate Z-scores for the current sequence against all sequences from Sample 1 and Sample 2
                z_scores_vs_d1 = []
                z_scores_vs_d2 = []

                for h_other in head_all:
                    if h_current == h_other: # Don't compare sequence to itself
                        continue

                    dist = hash_all.get(h_current, {}).get(h_other)
                    if dist is None:
                        # This might happen if hash_all wasn't fully populated, though unlikely with the loop structure
                        print(f"Warning: Missing distance between {h_current} and {h_other} during Z-score calculation.", file=sys.stderr)
                        continue

                    # Check which sample h_other belongs to
                    if h_other.startswith(s1_name + '_'):
                        z_scores_vs_d1.append((dist - mean_d1) / stdev_d1)
                    elif h_other.startswith(s2_name + '_'):
                        z_scores_vs_d2.append((dist - mean_d2) / stdev_d2)

                # Calculate mean Z-scores
                mean_z_vs_d1 = np.mean(z_scores_vs_d1) if z_scores_vs_d1 else float('inf') # Avg Z vs Sample 1 sequences
                mean_z_vs_d2 = np.mean(z_scores_vs_d2) if z_scores_vs_d2 else float('inf') # Avg Z vs Sample 2 sequences

                # Write Z-scores to stats file
                z1_str = f"{mean_z_vs_d1:.6f}" if np.isfinite(mean_z_vs_d1) else "Inf"
                z2_str = f"{mean_z_vs_d2:.6f}" if np.isfinite(mean_z_vs_d2) else "Inf"
                f_stats.write(f"{h_current}\t{z1_str}\t{z2_str}\n")

                # Determine classification based on the current sequence's origin
                if current_sample_name == s1_name:
                    # Classifying a sequence from Sample 1
                    q1 = 1 if mean_z_vs_d1 >= args.z_cutoff else 0 # High Z vs own sample? (Should be low) -> 1 = Different
                    q2 = 1 if mean_z_vs_d2 >= args.z_cutoff else 0 # High Z vs other sample? -> 1 = Different

                    if q1 == 0 and q2 == 0: # Low Z vs S1, Low Z vs S2 -> Belongs to S1, but also close to S2? (Mixed/AA)
                        results[s1_name]['AA'] += 1
                    elif q1 == 1 and q2 == 0: # High Z vs S1, Low Z vs S2 -> Doesn't belong to S1, close to S2? (BA - Belongs to S2?)
                        results[s1_name]['BA'] += 1
                    elif q1 == 0 and q2 == 1: # Low Z vs S1, High Z vs S2 -> Belongs to S1, different from S2 (AB - Homogeneous S1)
                        results[s1_name]['AB'] += 1
                    else: # q1 == 1 and q2 == 1 -> High Z vs S1, High Z vs S2 -> Doesn't belong to either? (BB - Unclassified)
                        results[s1_name]['BB'] += 1
                elif current_sample_name == s2_name:
                    # Classifying a sequence from Sample 2
                    # Note: Perl script compared S2 sequences vs S1 mean (q2_perl) and S2 sequences vs S2 mean (q1_perl)
                    q2_perl = 1 if mean_z_vs_d1 >= args.z_cutoff else 0 # High Z vs S1? -> 1 = Different from S1
                    q1_perl = 1 if mean_z_vs_d2 >= args.z_cutoff else 0 # High Z vs S2? (Should be low) -> 1 = Different from S2

                    if q1_perl == 0 and q2_perl == 0: # Low Z vs S2, Low Z vs S1 -> Belongs to S2, but also close to S1? (Mixed/AA)
                        results[s2_name]['AA'] += 1
                    elif q1_perl == 1 and q2_perl == 0: # High Z vs S2, Low Z vs S1 -> Doesn't belong to S2, close to S1? (AB - Belongs to S1?)
                        results[s2_name]['AB'] += 1
                    elif q1_perl == 0 and q2_perl == 1: # Low Z vs S2, High Z vs S1 -> Belongs to S2, different from S1 (BA - Homogeneous S2)
                        results[s2_name]['BA'] += 1
                    else: # q1_perl == 1 and q2_perl == 1 -> High Z vs S2, High Z vs S1 -> Doesn't belong to either? (BB - Unclassified)
                        results[s2_name]['BB'] += 1
            # End of classification loop inside 'with open'

    except IOError as e:
        print(f"Error writing stats file {stats_file_path}: {e}", file=sys.stderr)
    except Exception as e:
        # Catch other potential errors during classification/writing
        print(f"Error during classification or stats file writing: {e}", file=sys.stderr)
        # Ensure results dict exists even if loop failed partially
        if 'results' not in locals():
             results = {s1_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0},
                        s2_name: {'AA': 0, 'AB': 0, 'BA': 0, 'BB': 0}}

    # --- Determine Final Results ---
    print("Determining final classification...", file=sys.stderr)
    # (Classification logic remains the same as previous version)
    res1 = results[s1_name]

    # --- Determine Final Results ---
    print("Determining final classification...", file=sys.stderr)
    # (Classification logic remains the same as previous version)
    res1 = results[s1_name]
    tot1 = sum(res1.values())
    result_s1 = "Failed"
    if tot1 > 0:
        if res1['AA'] / tot1 >= 0.8:
            result_s1 = "Mixed"
        elif res1['AB'] / tot1 >= 0.5:
            result_s1 = "Homogenious"
        elif res1['BA'] / tot1 >= 0.5:
            result_s1 = "recipient"
        elif res1['BB'] > res1['AA'] and res1['BB'] > res1['AB'] and res1['BB'] > res1['BA']:
            result_s1 = "Too many unclassified reads"

    res2 = results[s2_name]
    tot2 = sum(res2.values())
    result_s2 = "Failed"
    if tot2 > 0:
        if res2['AA'] / tot2 >= 0.2:
            result_s2 = "Mixed"
        elif res2['AB'] / tot2 >= 0.5:
            result_s2 = "recipient"
        elif res2['BA'] / tot2 >= 0.5:
            result_s2 = "Homogenious"
        elif res2['BB'] > res2['AA'] and res2['BB'] > res2['AB'] and res2['BB'] > res2['BA']:
            result_s2 = "Too many unclassified reads"

    # --- Print Summary to STDOUT ---
    print("Sample\tIntermixture\tGroup1\tGroup2\tUnclustered")
    print("--------------------------------------------")
    print(f"{s1_name}:\t{res1['AA']}\t{res1['AB']}\t{res1['BA']}\t{res1['BB']}")
    print(f"{s2_name}:\t{res2['AA']}\t{res2['AB']}\t{res2['BA']}\t{res2['BB']}")
    print(f"\nIntermediate Classifications: {result_s1} / {result_s2}\n")

    # Determine final conclusion (Logic remains the same)
    conclusion = "Could Not Determine"
    if "Failed" in result_s1 or "Failed" in result_s2 or "unclassified" in result_s1 or "unclassified" in result_s2:
        conclusion = "Failed to confirm data transmission"
    elif result_s1 == "Mixed" and result_s2 == "Mixed":
        conclusion = "Complete intersample mixture"
    elif result_s1 == "Mixed" and result_s2 == "Homogenious":
        conclusion = f"Transmission from {s2_name} to {s1_name}"
    elif result_s2 == "Mixed" and result_s1 == "Homogenious":
        conclusion = f"Transmission from {s1_name} to {s2_name}"
    elif result_s1 == "Homogenious" and result_s2 == "Homogenious":
        conclusion = "No transmission detected"
    elif result_s1 == "recipient" and result_s2 == "Homogenious":
        conclusion = f"Transmission from {s2_name} to {s1_name}"
    elif result_s2 == "recipient" and result_s1 == "Homogenious":
        conclusion = f"Transmission from {s1_name} to {s2_name}"

    print(f"{s1_name} vs {s2_name}:")
    print(conclusion)
    print("______________________________________\n")

    # --- Generate R Plot ---
    print("Generating R plot...", file=sys.stderr)


    # 3. Generate R script only if legend generation succeeded
    # --- Generate t-SNE Plot using External R Script ---
    # This section assumes hash_all and head_all are populated from internal calculations
    print("Generating t-SNE plot using external R script...", file=sys.stderr)
    tsne_input_matrix_file = temp_dir / "ident.pmatrix" # Input filename for R script
    external_r_script_path = Path("/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline/plot_tsne.R") # Path to your external R script
    temp_pdf_output_name = "ident_tsne.pdf" # Expected output name from plot_tsne.R
    temp_pdf_output_path = temp_dir / temp_pdf_output_name
    final_pdf_output_path = reports_dir / f"Rtsne_{s1_name}_{s2_name}.pdf" # Final desired name and location
    print(f"Executing R script: {external_r_script_path} with input {tsne_input_matrix_file.name}", file=sys.stderr)
    # Use Rscript path from arguments
    r_ret = run_command([args.rscript_path, str(external_r_script_path), str(tsne_input_matrix_file)], cwd=temp_dir, log_file=log_file)
    
    # 4. Check result and move/rename the PDF
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


if __name__ == "__main__":
        # 1. Run snp-dists to generate distance matrix
    pmatrix_file = temp_dir / "ident.pmatrix"
    # Construct the full command string with shell redirection
    # Ensure paths with spaces are handled if necessary (though unlikely here)
    # Using f-string for clarity, but ensure args.snp_dists_path doesn't need complex quoting
    command_str = f"{args.snp_dists_path} -m {data3_unique_headers_out} > {pmatrix_file}"
    
    print(f"Running command via shell: {command_str}", file=sys.stderr)
    subprocess.run(command_str, shell=True, cwd=temp_dir)
    main()