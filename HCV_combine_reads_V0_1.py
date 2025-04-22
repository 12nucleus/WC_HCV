# HCV_combine_reads_V0_1.py
import argparse
import subprocess
import gzip
from pathlib import Path
from collections import Counter
import sys
import re
import os
import shutil

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def run_command(cmd_list, log_file=None, error_log_file=None, cwd=None):
    """Runs an external command and handles output."""
    print(f"Running command: {' '.join(cmd_list)}", file=sys.stderr)
    stdout_pipe = subprocess.PIPE
    stderr_pipe = subprocess.PIPE
    if log_file:
        stdout_pipe = open(log_file, 'w')
    if error_log_file:
        stderr_pipe = open(error_log_file, 'w')

    try:
        process = subprocess.run(cmd_list, check=True, text=True,
                                 stdout=stdout_pipe, stderr=stderr_pipe, cwd=cwd)
        if log_file:
            stdout_pipe.close()
        if error_log_file:
            stderr_pipe.close()
        return process
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd_list)}", file=sys.stderr)
        print(f"Return code: {e.returncode}", file=sys.stderr)
        if not log_file:
             print(f"Stdout: {e.stdout}", file=sys.stderr)
        if not error_log_file:
             print(f"Stderr: {e.stderr}", file=sys.stderr)
        if log_file:
            stdout_pipe.close()
        if error_log_file:
            stderr_pipe.close()
        sys.exit(1) # Exit if a critical command fails
    except FileNotFoundError:
        print(f"Error: Command not found: {cmd_list[0]}. Ensure it's installed and in PATH.", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Combine paired-end reads, filter based on k-mers, and generate k-mer counts.")
    parser.add_argument("working_dir", type=Path, help="Working directory containing input and for output subdirectories.")
    parser.add_argument("r1_file", type=Path, help="Path to the R1 FASTQ file.")
    parser.add_argument("tmp_dir", type=Path, help="Path to the shared temporary directory for this run.")
    parser.add_argument("--ref_kmer_file", type=Path, default=Path("/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline/reference_genome/"),
                        help="Path to the reference k-mer file (default: /Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline/reference_genome/).")
    parser.add_argument("--cutadapt_path", type=str, default="cutadapt", help="Path to the cutadapt executable.")
    parser.add_argument("--bbmerge_path", type=str, default="/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline/bbmap/bbmerge.sh", help="Path to the bbmerge.sh script.")
    parser.add_argument("--min_overlap", type=int, default=100, help="Minimum overlap for bbmerge.")
    parser.add_argument("--min_count", type=int, default=15, help="Minimum count for a sequence to be kept.")
    parser.add_argument("--min_kmer_matches", type=int, default=15, help="Minimum number of reference k-mer matches required.")
    parser.add_argument("--kmer_size_jf", type=int, default=125, help="K-mer size for jellyfish counting.")
    parser.add_argument("--min_final_reads", type=int, default=50, help="Minimum number of filtered reads required to PASS.")
    parser.add_argument("--keep_tmp", action='store_true', help="Keep temporary files.")
    parser.add_argument("--keep_unmerged", action='store_true', help="Keep the unmerged R1 and R2 reads from bbmerge.")

    args = parser.parse_args()

    # --- Validate inputs and create directories ---
    if not args.working_dir.is_dir():
        print(f"Error: Working directory not found: {args.working_dir}", file=sys.stderr)
        sys.exit(1)
    if not args.r1_file.is_file():
        print(f"Error: R1 file not found: {args.r1_file}", file=sys.stderr)
        sys.exit(1)
    if not args.tmp_dir.is_dir():
        print(f"Error: Temporary directory not found: {args.tmp_dir}", file=sys.stderr)
        sys.exit(1)

    combined_reads_dir = args.working_dir / "combined_reads"
    fasta_dir = args.working_dir / "HCV_fasta"
    kmers_dir = args.working_dir / "HCV_Kmers"
    combined_reads_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)
    kmers_dir.mkdir(parents=True, exist_ok=True)

    # Validate the reference k-mer file/directory using hardcoded values
    ref_kmer_dir = Path("/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline/reference_genome/")
    ref_file1 = ref_kmer_dir / "1a_kmer_ref_counts.out"
    ref_file2 = ref_kmer_dir / "1b_kmer_ref_counts.out"
    if not (ref_file1.is_file() and ref_file2.is_file()):
        print(f"Error: Expected reference k-mer files not found in {ref_kmer_dir}", file=sys.stderr)
        sys.exit(1)
    ref_files = [ref_file1, ref_file2]

    kmer_check = set()
    try:
        for ref_file in ref_files:
            with open(ref_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        continue
                    kmer = line.strip().upper()
                    if kmer:
                        kmer_check.add(kmer)
    except IOError as e:
        print(f"Error reading reference k-mer file(s): {e}", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(kmer_check)} reference k-mers.", file=sys.stderr)

    # --- Determine R2 filename and prefix ---
    r1_filename = args.r1_file.name
    try:
        # Attempt common naming conventions
        if "_R1_" in r1_filename:
            r2_filename = r1_filename.replace("_R1_", "_R2_", 1)
        elif "_R1." in r1_filename:
             r2_filename = r1_filename.replace("_R1.", "_R2.", 1)
        elif "_1." in r1_filename:
             r2_filename = r1_filename.replace("_1.", "_2.", 1)
        elif "_1_" in r1_filename:
             r2_filename = r1_filename.replace("_1_", "_2_", 1)
        else:
            raise ValueError("Cannot determine R2 filename pattern from R1 filename.")

        r2_file = args.r1_file.parent / r2_filename
        if not r2_file.is_file():
             print(f"Error: Deduced R2 file not found: {r2_file}", file=sys.stderr)
             sys.exit(1)

        # Extract prefix (assuming format like PREFIX_...)
        prefix = r1_filename.split('_')[0]
        if not prefix:
             prefix = args.r1_file.stem # Fallback if no underscore

    except Exception as e:
        print(f"Error determining R2 filename or prefix: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing sample: {prefix}", file=sys.stderr)
    print(f"R1 file: {args.r1_file}", file=sys.stderr)
    print(f"R2 file: {r2_file}", file=sys.stderr)

    # --- Define intermediate and final file paths ---
    # Intermediate files go into tmp_dir
    trimmed_r1 = args.tmp_dir / f"{prefix}.1.fastq"
    trimmed_r2 = args.tmp_dir / f"{prefix}.2.fastq"
    jf_output = args.tmp_dir / f"{prefix}_mer_counts.jf" # Jellyfish intermediate file
    cutadapt_log = args.tmp_dir / f"{prefix}_cutadapt.log"
    bbmerge_log = args.tmp_dir / f"{prefix}_bbmerge.log"
    # Jellyfish dump output (before moving) also goes to tmp_dir
    kmer_dump_output = args.tmp_dir / f"{prefix}.kmers"

    # Final outputs go to designated directories relative to working_dir
    merged_fq_gz = combined_reads_dir / f"{prefix}.fq.gz" # BBmerge output is final
    final_fasta = fasta_dir / f"{prefix}.fasta"
    final_fasta_gz = fasta_dir / f"{prefix}.fasta.gz"
    # Final kmer path (after moving)
    final_kmer_output_gz = kmers_dir / f"{prefix}.kmers.gz"
    unmerged_r1_fq_gz = combined_reads_dir / f"{prefix}.unmerged.1.fq.gz"
    unmerged_r2_fq_gz = combined_reads_dir / f"{prefix}.unmerged.2.fq.gz"


    # --- Run Cutadapt ---
    # Original primers: -G GGATATGATGATGAACTGGT -g ATGTGCCAGCTGCCGTTGGTGT -g GGATATGATGATGAACTGGT -G ATGTGCCAGCTGCCGTTGGTGT
    # Assuming these are 5' adapters on R1/R2 respectively, and maybe also 3'? Let's simplify based on common use.
    # Using -a for R1 3' adapter and -A for R2 3' adapter. The original used -g/-G which are 5'. Re-check if this is correct.
    # The original command seems to specify the *same* adapters twice with -g and -G. This might be removing 5' adapters.
    # Let's stick to the original -g/-G usage.
    cutadapt_cmd = [
        args.cutadapt_path,
        "-G", "GGATATGATGATGAACTGGT", # 5' adapter Read 1
        "-g", "ATGTGCCAGCTGCCGTTGGTGT", # 5' adapter Read 2? Seems reversed based on common pairing. Let's assume original script knew best.
        "-g", "GGATATGATGATGAACTGGT", # 5' adapter Read 2 again?
        "-G", "ATGTGCCAGCTGCCGTTGGTGT", # 5' adapter Read 1 again? This is confusing. Sticking to original.
        "-o", str(trimmed_r1),
        "-p", str(trimmed_r2),
        str(args.r1_file),
        str(r2_file)
    ]
    print("Running Cutadapt...", file=sys.stderr)
    run_command(cutadapt_cmd, log_file=cutadapt_log)

    # --- Run BBmerge ---
    bbmerge_cmd = [
        args.bbmerge_path,
        f"pfilter=1",
        f"minoverlap0={args.min_overlap}", # Allow no mismatches in initial overlap seed
        f"minoverlap={args.min_overlap}", # Minimum overlap length
        f"in1={trimmed_r1}",
        f"in2={trimmed_r2}",
        f"out={merged_fq_gz}",
        f"ziplevel=6"
    ]
    if args.keep_unmerged:
        bbmerge_cmd.extend([
            f"outu1={unmerged_r1_fq_gz}",
            f"outu2={unmerged_r2_fq_gz}"
        ])
        print("BBmerge will keep unmerged reads.", file=sys.stderr)

    # stderr is captured for parsing
    print("Running BBmerge...", file=sys.stderr)
    bbmerge_process = run_command(bbmerge_cmd, error_log_file=bbmerge_log) # BBmerge writes stats to stderr

    # --- Parse BBmerge Log ---
    total_pairs_bb = "N/A"
    percent_merged_bb = "N/A"
    try:
        with open(bbmerge_log, 'r') as log_f:
            for line in log_f:
                if 'Pairs:' in line:
                    parts = line.strip().split('\t')
                    if len(parts) > 1:
                        total_pairs_bb = parts[1]
                elif 'Joined:' in line:
                     parts = line.strip().split('\t')
                     if len(parts) > 2:
                        percent_merged_bb = parts[2]
    except IOError as e:
        print(f"Warning: Could not read bbmerge log {bbmerge_log}: {e}", file=sys.stderr)
    except Exception as e:
         print(f"Warning: Could not parse bbmerge log {bbmerge_log}: {e}", file=sys.stderr)

    print(f"Total reads\tPercent merged\t# HVR1\t# Rejected\tStatus") # Header
    print(f"{total_pairs_bb}\t{percent_merged_bb}\t", end="") # Print initial stats

    # --- Cleanup Intermediate Files (Trimmed Reads) ---
    # Logs and jf files are cleaned up later if PASS and not keep_tmp
    if not args.keep_tmp:
        print("Cleaning up trimmed read files...", file=sys.stderr)
        try:
            trimmed_r1.unlink(missing_ok=True)
            trimmed_r2.unlink(missing_ok=True)
        except OSError as e:
             print(f"Warning: Could not delete trimmed read files: {e}", file=sys.stderr)


    # --- Process Merged Reads ---
    print(f"Processing merged reads from {merged_fq_gz}...", file=sys.stderr)
    seq_counts = Counter()
    total_reads_processed = 0
    try:
        with gzip.open(merged_fq_gz, "rt") as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip().upper()
                plus = f.readline()
                qual = f.readline()
                if not header.startswith('@') or not plus.startswith('+') or not seq:
                     print(f"Warning: Malformed FASTQ entry near read {total_reads_processed+1} in {merged_fq_gz}", file=sys.stderr)
                     continue # Skip malformed entry

                total_reads_processed += 1
                # No length filter applied here as per the commented-out Perl code
                # No primer trimming applied here as per the commented-out Perl code
                seq_counts[seq] += 1
    except FileNotFoundError:
         print(f"Error: Merged file not found: {merged_fq_gz}", file=sys.stderr)
         sys.exit(1)
    except Exception as e:
        print(f"Error processing merged FASTQ file {merged_fq_gz}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(seq_counts)} unique sequences from {total_reads_processed} merged reads.", file=sys.stderr)

    # --- Filter Sequences ---
    print(f"Filtering sequences (min_count={args.min_count}, min_kmer_matches={args.min_kmer_matches})...", file=sys.stderr)
    filtered_sequences = []
    total_filtered_reads = 0
    rejected_min_representation = 0
    rejected_kmer_matches = 0
    tag = 1

    # Sort by count descending (optional, matches Perl script)
    #sorted_seqs = sorted(seq_counts.items(), key=lambda item, count: item[1], reverse=True)
    sorted_seqs = sorted(seq_counts.items(), key=lambda item: item[1], reverse=True)
    for seq, count in sorted_seqs:
        if count >= args.min_count:
            match_count = 0
            seq_rc = None # Calculate only if needed
            for ref_kmer in kmer_check:
                if ref_kmer in seq:
                    match_count += 1
                    continue # Found forward match
                # Check reverse complement only if forward not found
                if seq_rc is None:
                    seq_rc = reverse_complement(seq)
                if ref_kmer in seq_rc:
                    match_count += 1

            if match_count >= args.min_kmer_matches:
                total_filtered_reads += count
                header = f">{tag}_{prefix}_{count}"
                filtered_sequences.append((header, seq))
            else:
                rejected_kmer_matches += count
        else:
            rejected_min_representation += count
        tag += 1

    print(f"{total_filtered_reads}\t{rejected_kmer_matches}\t", end="") # Print filtering stats

    # --- Final Output and Status ---
    if total_filtered_reads >= args.min_final_reads:
        print("PASS") # Print status to stdout
        print(f"Writing {len(filtered_sequences)} sequences to {final_fasta}...", file=sys.stderr)
        try:
            with open(final_fasta, 'w') as f_out:
                for header, seq in filtered_sequences:
                    f_out.write(f"{header}\n{seq}\n")
        except IOError as e:
             print(f"Error writing FASTA file {final_fasta}: {e}", file=sys.stderr)
             # Continue to jellyfish/gzip if possible, but report error

        # --- Run Jellyfish in Docker ---
        # Mount working directory and output directory as volumes
        docker_cmd = [
            "docker", "run", "--rm", # Remove container after run
            "-v", f"{fasta_dir.resolve()}:/data:ro", # Mount final fasta dir read-only
            "-v", f"{args.tmp_dir.resolve()}:/output", # Mount tmp_dir for output
            "dnalinux/jellyfish",
            "jellyfish", "count",
            "-m", str(args.kmer_size_jf),
            "-s", "100M", # Assuming this is memory allocation
            "-t", "1",    # Assuming single thread
            "-o", "/output/mer_counts.jf", # Output inside container
            "/data/" + final_fasta.name # Input inside container
        ]
        print("Running Jellyfish count in Docker...", file=sys.stderr)
        run_command(docker_cmd)
        docker_dump_cmd = [
            "docker", "run", "--rm",
            "-v", f"{args.tmp_dir.resolve()}:/output", # Mount tmp_dir
            "dnalinux/jellyfish",
            "jellyfish", "dump",
            "/output/mer_counts.jf", # Input jf file inside container
            "-o", f"/output/{kmer_dump_output.name}" # Output dump file inside container (in tmp_dir)
        ]
        print("Running Jellyfish dump in Docker...", file=sys.stderr)
        run_command(docker_dump_cmd)

        # --- Gzip Outputs ---
        print("Compressing output files...", file=sys.stderr)
        gzip_fasta_cmd = ["gzip", "--force", str(final_fasta)]
        run_command(gzip_fasta_cmd)

        # Correct path: Jellyfish dump output is in args.tmp_dir
        # kmer_file_to_compress = args.output_dir / kmer_output.name # Old logic
        gzip_kmers_cmd = ["gzip", "--force", str(kmer_dump_output)] # Compress the file in tmp_dir
        run_command(gzip_kmers_cmd)

        # Move the compressed kmer file from tmp_dir to the final destination
        compressed_kmer_source = args.tmp_dir / f"{kmer_dump_output.name}.gz" # e.g., tmp_dir/prefix.kmers.gz
        final_kmer_dest = final_kmer_output_gz # e.g., HCV_Kmers/prefix.kmers.gz
        try:
            print(f"Moving {compressed_kmer_source} to {final_kmer_dest}", file=sys.stderr)
            shutil.move(str(compressed_kmer_source), str(final_kmer_dest))
        except Exception as e:
            print(f"Error moving kmer file {compressed_kmer_source} to {final_kmer_dest}: {e}", file=sys.stderr)
            # Decide if this should cause a failure or just a warning
            # sys.exit(1) # Uncomment to make this error fatal

        # --- Final Cleanup (if PASS and not keep_tmp) ---
        if not args.keep_tmp:
            print("Cleaning up intermediate files from tmp dir...", file=sys.stderr)
            files_to_delete = [jf_output, bbmerge_log, cutadapt_log, kmer_dump_output] # kmer_dump_output is the uncompressed one
            for f_path in files_to_delete:
                try:
                    f_path.unlink(missing_ok=True)
                except OSError as e:
                    print(f"Warning: Could not delete temporary file {f_path}: {e}", file=sys.stderr)

    else:
        print("FAIL") # Print status to stdout
        print(f"FAIL: Total filtered reads ({total_filtered_reads}) below threshold ({args.min_final_reads}).", file=sys.stderr)

if __name__ == "__main__":
    main()
