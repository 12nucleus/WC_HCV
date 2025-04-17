# HCV_kmers_distancer_V0_1.py
import argparse
import gzip
from pathlib import Path
import sys
import os

def load_kmers(kmer_file_path: Path) -> set:
    """Loads unique k-mers from a gzipped file (like jellyfish dump output)."""
    if not kmer_file_path.is_file():
        raise FileNotFoundError(f"K-mer file not found: {kmer_file_path}")

    kmers = set()
    try:
        with gzip.open(kmer_file_path, "rt") as f:
            for line in f:
                # Jellyfish dump format often includes >header lines
                # This script assumes simple list of kmers or ignores headers
                if not line.startswith('>'):
                    kmer = line.strip().upper()
                    if kmer: # Avoid adding empty lines if present
                        kmers.add(kmer)
    except Exception as e:
        print(f"Error reading k-mer file {kmer_file_path}: {e}", file=sys.stderr)
        # Return empty set or re-raise depending on desired behavior
        # For this script, continuing with an empty set might be okay,
        # but it's better to signal the error.
        raise # Re-raise the exception
    return kmers

def main():
    parser = argparse.ArgumentParser(description="Calculate k-mer overlap between a query sample and other samples.")
    parser.add_argument("working_dir", type=Path, help="Working directory containing the HCV_Kmers subdirectory.")
    parser.add_argument("query_prefix", type=str, help="Prefix of the query k-mer file (e.g., 'sampleA' for 'sampleA.kmers.gz').")
    parser.add_argument("--kmer_dir_name", type=str, default="HCV_Kmers", help="Name of the subdirectory containing kmer files.")
    parser.add_argument("--overlap_threshold", type=float, default=0.05, help="Threshold for marking significant overlap.")

    args = parser.parse_args()

    kmer_dir = args.working_dir / args.kmer_dir_name
    if not kmer_dir.is_dir():
        print(f"Error: K-mer directory not found: {kmer_dir}", file=sys.stderr)
        sys.exit(1)

    query_kmer_file = kmer_dir / f"{args.query_prefix}.kmers.gz"

    # --- Load Query K-mers ---
    try:
        print(f"Loading query k-mers from: {query_kmer_file}", file=sys.stderr)
        query_kmers = load_kmers(query_kmer_file)
        total_query = len(query_kmers)
        print(f"Loaded {total_query} unique k-mers for query {args.query_prefix}.", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Query k-mer file not found: {query_kmer_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading query k-mers: {e}", file=sys.stderr)
        sys.exit(1)


    # --- Print Header ---
    print("S1\tS2\tKmer Overlaps\tTotal S1\tTotal S2\tRatio")

    # --- Iterate Through Subject K-mer Files ---
    subject_files = list(kmer_dir.glob("*.kmers.gz"))
    if not subject_files:
        print(f"Warning: No '.kmers.gz' files found in {kmer_dir}", file=sys.stderr)

    for subject_file in subject_files:
        # Extract prefix (robustly handles potential dots in prefix)
        subject_prefix = subject_file.name.replace(".kmers.gz", "")

        if subject_prefix == args.query_prefix:
            continue # Skip self-comparison

        try:
            # --- Load Subject K-mers ---
            # print(f"Loading subject k-mers from: {subject_file}", file=sys.stderr) # Verbose
            subject_kmers = load_kmers(subject_file)
            total_subject = len(subject_kmers)
            # print(f"Loaded {total_subject} unique k-mers for subject {subject_prefix}.", file=sys.stderr) # Verbose

            # --- Calculate Overlap ---
            overlap_kmers = query_kmers.intersection(subject_kmers)
            overlap_count = len(overlap_kmers)

            # --- Calculate Ratio and Print ---
            if overlap_count >= 1:
                # Ensure no division by zero if a file somehow had 0 unique kmers
                min_total = min(total_query, total_subject)
                if min_total > 0:
                    ratio = overlap_count / min_total
                else:
                    ratio = 0.0 # Define ratio as 0 if one set is empty

                marker = "*" if ratio >= args.overlap_threshold else ""
                print(f"{args.query_prefix}\t{subject_prefix}\t{overlap_count}\t{total_query}\t{total_subject}\t{ratio:.5f}{marker}")

        except FileNotFoundError:
            # Should not happen if glob worked, but good practice
            print(f"Warning: Subject file disappeared?: {subject_file}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"Warning: Error processing subject file {subject_file}: {e}", file=sys.stderr)
            continue # Skip this subject file on error

    # --- Print Footer ---
    print(f"* indicated Kmer overlap ratio >= {args.overlap_threshold}")

if __name__ == "__main__":
    main()