# HCV Transmission Analysis Pipeline (Python Conversion)

This pipeline analyzes Hepatitis C Virus (HCV) sequencing data to identify potential transmission links between samples. It is a Python conversion of an original set of Perl scripts.

## Pipeline Overview

The pipeline performs the following main steps:

1.  **Read Processing (`HCV_combine_reads_V0_1.py`):**
    *   Trims adapters from paired-end FASTQ reads using `cutadapt`.
    *   Merges overlapping reads using `bbmerge.sh`.
    *   Filters merged reads based on minimum count and presence of reference HVR1 k-mers.
    *   Generates a FASTA file of filtered sequences and calculates k-mer counts using `jellyfish` (run inside a Docker container).
2.  **K-mer Comparison (`HCV_kmers_distancer_V0_1.py`):**
    *   Compares the k-mer sets generated for a query sample against all other samples.
    *   Calculates an overlap ratio to identify samples with significant k-mer sharing (potential links).
3.  **Transmission Testing (`HCV_transmission_test_V0_1.py`):**
    *   Takes pairs of samples identified as potential links by the k-mer comparison.
    *   Aligns the filtered FASTA sequences within and between the pair using `pegi3s/mafft` Docker image.
    *   Calculates pairwise sequence distances using a custom metric (optimized with NumPy).
    *   Performs Z-score analysis based on within-group and between-group distances to classify the relationship (e.g., transmission, mixture, unrelated).
    *   Generates a distance matrix using `snp-dists`.
    *   Creates a t-SNE visualization plot using R (`Rtsne`, `ggplot2`).
4.  **Orchestration (`HCV_Master_script_V0_1.py`):**
    *   Coordinates the execution of the above steps for all input samples.
    *   Manages input/output files and temporary directories.
    *   Generates summary reports for each sample.

## Code Layout

The pipeline consists of the following Python scripts:

*   `HCV_Master_script_V0_1.py`: The main script to run the entire pipeline.
*   `HCV_combine_reads_V0_1.py`: Handles adapter trimming, read merging, filtering, and k-mer generation.
*   `HCV_kmers_distancer_V0_1.py`: Compares k-mer sets between samples.
*   `HCV_transmission_test_V0_1.py`: Performs detailed transmission analysis using alignments and generates plots.

## Dependencies

### 1. Python Environment
*   Python 3.x
*   A Conda environment (e.g., named `HCV`) is recommended for managing dependencies.

### 2. Python Packages
Install these within your activated conda environment using pip:
```bash
pip install numpy biopython
```

### 3. Bioinformatics Tools
These tools must be installed and accessible in your system's `PATH`. **Do not install these using conda** as per user instructions.
*   `cutadapt`: For adapter trimming.
*   `bbmap`: Provides `bbmerge.sh` for read merging.
*   `snp-dists`: For generating SNP distance matrices from alignments. brew install brewsci/bio/snp-dists

### 4. Docker
*   [Docker](https://www.docker.com/) must be installed and running.

### 5. R Environment
*   R installation accessible from the command line.
*   R packages: `Rtsne`, `ggplot2`.

## Setup

1.  **Create and Activate Conda Environment (Recommended):**
    ```bash
    conda create -n HCV python=3.9 # Or your preferred Python 3 version
    conda activate HCV
    ```

2.  **Install Python Packages:**
    ```bash
    pip install numpy biopython
    ```

3.  **Install Bioinformatics Tools:**
    Ensure `cutadapt`, `bbmerge.sh` (from bbmap), and `snp-dists` are installed system-wide or otherwise available in your `PATH`. Verify their installation (e.g., `cutadapt --version`, `snp-dists --version`).

4.  **Install Docker:**
    [Docker](https://www.docker.com/) must be installed and running to use `jellyfish` and `mafft`.

5.  **Install R Packages:**
    *   Start R:
        ```bash
        R
        ```
    *   Inside the R console, run:
        ```R
        install.packages(c("Rtsne", "ggplot2"), repos='http://cran.us.r-project.org')
        ```
    *   Quit R:
        ```R
        q()
        # Choose 'n' if asked to save workspace image
        ```

## Running the Pipeline

The pipeline is executed using the master script `HCV_Master_script_V0_1.py`.

**Basic Usage:**

```bash
conda activate HCV
python HCV_Master_script_V0_1.py <path_to_reads_directory> <path_to_pipeline_base_directory>
```

**Arguments:**

*   `<path_to_reads_directory>`: (Required) Path to the directory containing the input FASTQ files (e.g., `*_R1_001.fastq.gz`).
*   `<path_to_pipeline_base_directory>`: (Required) Path to the directory where the Python scripts (`.py`) are located. Output subdirectories (`Reports`, `HCV_Kmers`, `HCV_fasta`, `combined_reads`) will be created here.

**Example:**

If your scripts are in `/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline` and your reads are in `/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline/reads`:

```bash
cd /Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline
conda activate HCV
python HCV_Master_script_V0_1.py ./reads .
```

**Optional Arguments (`HCV_Master_script_V0_1.py`):**

Run `python HCV_Master_script_V0_1.py -h` to see all options, including:

*   `--keep_tmp`: Keep temporary directories created during runs (useful for debugging).
*   `--combine_script`, `--distancer_script`, `--transmission_script`: Specify alternative names for the component scripts if needed.
*   `--r1_pattern`: Change the glob pattern used to find R1 FASTQ files (default: `*_R1_001.fastq.gz`).
*   `--kmer_overlap_threshold`: Adjust the k-mer overlap ratio threshold for triggering the transmission test (default: 0.05).

The `jellyfish` k-mer counting and `mafft` alignment steps are run inside Docker containers. Ensure Docker is installed and running.

**Optional Arguments (Component Scripts):**

The master script passes necessary arguments to the component scripts. However, the component scripts also have their own optional arguments (e.g., paths to executables like `snp-dists` if not in PATH, specific parameters like minimum overlap). If you need to modify these, you would typically adjust the defaults within the respective Python script files or modify the `HCV_Master_script_V0_1.py` to pass them through.

## Output

*   **`Reports/`**: Contains text reports (`<sample>_report.out`) summarizing the processing steps and results for each sample. Also contains PDF plots (`Rtsne_<sample1>_<sample2>.pdf`) visualizing the t-SNE analysis for potentially linked pairs.
*   **`combined_reads/`**: Contains merged FASTQ files (`<sample>.fq.gz`) generated by `bbmerge`.
*   **`HCV_fasta/`**: Contains gzipped FASTA files (`<sample>.fasta.gz`) of the filtered sequences passing the `HCV_combine_reads` step.
*   **`HCV_Kmers/`**: Contains gzipped k-mer count files (`<sample>.kmers.gz`) generated by `jellyfish`.
*   **Standard Output/Error:** Progress messages, classifications, and potential warnings/errors are printed to the console during execution.