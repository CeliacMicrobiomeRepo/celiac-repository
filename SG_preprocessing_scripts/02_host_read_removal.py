"""
Performs host read removal from paired-end shotgun metagenomic sequencing data
using Bowtie2.

This script iterates through specified dataset directories, identifies paired-end
FASTQ files, aligns them against a host reference genome index using Bowtie2,
and outputs the unaligned (non-host) reads to a designated directory.
It also generates a log file summarizing the process and the percentage of
reads removed for each sample. It also creates a TSV summary file within each
dataset directory detailing the removal percentage per sample.

Assumes:
- Input FASTQ files are paired-end and reside in a 'input_fastqs' subdirectory
  within each dataset directory.
- FASTQ files are gzipped (.fq.gz).
- A consistent naming convention for paired reads (e.g., containing '_1' and '_2').
- Bowtie2 executable is in the system PATH or its full path is provided.
- A pre-built Bowtie2 index for the host genome exists.

Directory Structure:
Root/
├── SG_preprocessing_scripts/
│   └── 02_host_read_removal.py         <-- This script
├── SG_datasets/
│   ├── <Dataset_ID_1>/                 # e.g., SG_80_Mouzan
│   │   ├── input_fastqs/                     # Input trimmed FASTQs (paired)
│   │   │   ├── sampleA_1.fq.gz
│   │   │   ├── sampleA_2.fq.gz
│   │   │   └── ...
│   │   ├── host_removed_fastqs/        # Output directory (script creates this)
│   │   │   ├── sampleA_1.fq.gz
│   │   │   ├── sampleA_2.fq.gz
│   │   │   ├── sampleA_2_bowtie2.log
│   │   │   └── ...
│   │   └── about_host_read_removal.tsv      # Output TSV summary (script creates this)
│   └── <Dataset_ID_2>/
│       ├── input_fastqs/
│       │   └── ...
│       ├── host_removed_fastqs/
│       │   └── ...
│       └── about_host_read_removal.tsv
├── host_genome/                        # Directory for the host genome index
│   ├── GRCh38_noalt_decoy_as.1.bt2
│   ├── GRCh38_noalt_decoy_as.2.bt2
│   ├── ... (all bowtie2 index files)
└── host_removal_summary.log            # Overall log file (script creates this)


Set Up Bowtie2 and Python in Conda Environment (ubuntu):
Create a new conda environment for metagenomic analysis:
```
conda create -n metagenomics python=3.9
conda activate metagenomics
```
Install Bowtie2 and other dependencies:
```
conda install -c bioconda bowtie2 samtools
pip install biopython
```
Verify installation:
```
bowtie2 --version
python --version
```

Run the script (WSL):
```
conda activate metagenomics
cd /mnt/b/CMR_v1/SG_preprocessing_scripts
python 02_host_read_removal.py
```

Get Host Genome Index:
- Get the index files from the Bowtie2 website (https://benlangmead.github.io/aws-indexes/bowtie)
- Download the Human GRCh38 no-alt +decoy set1 (multiple `.bt2` files).
- Place them in your `host_genome` directory.
- Ensure the `HOST_INDEX_PATH` points to the base name (e.g., `C:\\path\\to\\host_genome\\GRCh38_noalt_decoy_as`, not including `.1.bt2`).


"""

import os
import subprocess
import glob
import logging
import re
import sys
from pathlib import Path

# ================== CONFIGURATION ==================

# --- Paths ---
# List of dataset directories to process
# Use absolute paths or paths relative to the script location
BASE_PROJECT_DIR = Path(__file__).resolve().parents[1] # Assumes script is in SG_preprocessing_scripts
DATASET_BASE_DIR = BASE_PROJECT_DIR / "SG_datasets"
# Specify dataset IDs relative to DATASET_BASE_DIR, or leave empty to process all dirs found
# Example: DATASET_IDS = ["SG_80_Mouzan", "SG_132_Francavilla"]
DATASET_IDS = ["SG_4_Example"] # ["SG_80_Mouzan", "SG_132_Francavilla"]

# Subdirectory within each dataset dir containing input FASTQs
INPUT_FASTQS_SUBDIR = "input_fastqs"
# Subdirectory within each dataset dir for output non-host FASTQs
OUT_FASTQS_SUBDIR = "host_removed_fastqs"
# Log file path
LOG_FILE = BASE_PROJECT_DIR / "host_removal_summary.log"

# --- Host Genome ---
# Path to the Bowtie2 index files for the host genome (excluding trailing extensions like .1.bt2)
# Example: HOST_INDEX_PATH = "/path/to/host_genome/GRCh38_noalt_decoy_as"
HOST_INDEX_PATH = BASE_PROJECT_DIR / "host_genome/GRCh38_noalt_decoy_as" # Adjust as needed

# --- Tools ---
# Full path to bowtie2 executable if not in PATH, otherwise just 'bowtie2'
BOWTIE2_EXE = "bowtie2"

# --- Parameters ---
# Number of threads for Bowtie2
NUM_THREADS = 26
# Bowtie2 alignment sensitivity preset (e.g., --very-fast-local, --sensitive-local, --very-sensitive-local)
BOWTIE2_SENSITIVITY = "--sensitive-local"
# Naming convention patterns to identify read pairs
# Adjust if your files use different patterns (e.g., _1.fq.gz, _2.fq.gz)
R1_PATTERN = "_1"
R2_PATTERN = "_2"

# ================== END CONFIGURATION ==================


# --- Logging Setup ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE, mode='w'),
        logging.StreamHandler(sys.stdout) # Also print logs to console
    ]
)

# --- Helper Functions ---
def find_fastq_pairs(fastq_dir, r1_pattern, r2_pattern):
    """Finds pairs of FASTQ files based on naming patterns."""
    pairs = {}
    r1_files = glob.glob(os.path.join(fastq_dir, f"*{r1_pattern}*.fq.gz"))
    logging.info(f"Found {len(r1_files)} potential R1 files in {fastq_dir}")

    for r1_file in r1_files:
        # Derive sample name and expected R2 file name
        base_name = os.path.basename(r1_file)
        try:
            sample_name = base_name.split(r1_pattern)[0]
            r2_file_name = base_name.replace(r1_pattern, r2_pattern, 1)
            r2_file = os.path.join(fastq_dir, r2_file_name)
        except Exception as e:
             logging.warning(f"Could not determine sample name or R2 filename for {r1_file}: {e}. Skipping.")
             continue

        if os.path.exists(r2_file):
            if sample_name in pairs:
                logging.warning(f"Duplicate sample name detected: {sample_name}. Overwriting previous entry.")
            pairs[sample_name] = (r1_file, r2_file)
            logging.debug(f"Paired {sample_name}: R1={r1_file}, R2={r2_file}")
        else:
            logging.warning(f"Could not find matching R2 file for {r1_file} (expected: {r2_file_name}). Skipping pair.")

    logging.info(f"Identified {len(pairs)} valid FASTQ pairs.")
    return pairs

def run_bowtie2(sample_name, r1_in, r2_in, host_index, out_dir, threads, sensitivity, bowtie2_exe):
    """Runs Bowtie2 to remove host reads."""
    output_prefix = os.path.join(out_dir, sample_name)
    unaligned_r1_out = f"{output_prefix}{R1_PATTERN}.fq.gz"
    unaligned_r2_out = f"{output_prefix}{R2_PATTERN}.fq.gz"
    log_path = f"{output_prefix}_bowtie2.log" # Log file for bowtie2 stderr

    # Using --un-conc-gz directly writes unaligned reads to FASTQ files
    cmd = [
        bowtie2_exe,
        "-x", str(host_index),
        "-1", r1_in,
        "-2", r2_in,
        "-S", "/dev/null", # Discard SAM alignment output
        "--threads", str(threads),
        sensitivity,
        "--no-unal", # Don't report unaligned reads in SAM (saves computation)
        "--un-conc-gz", output_prefix + "_%.fq.gz" # Template for unaligned output pairs
    ]
    # Note: The actual output files will be named based on the template,
    # e.g., sample_name_1.fq.gz and sample_name_2.fq.gz.
    # We will rename them later to match R1/R2 pattern.

    logging.info(f"Running Bowtie2 for {sample_name}...")
    logging.debug(f"Command: {' '.join(cmd)}")

    try:
        # Capture stderr for parsing alignment stats
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"Bowtie2 completed successfully for {sample_name}.")
        # Write Bowtie2's stderr (summary stats) to a log file
        with open(log_path, 'w') as log_f:
            log_f.write(result.stderr)
        logging.info(f"Bowtie2 log saved to: {log_path}")

        # --- Rename output files ---
        # Bowtie2 outputs sample_name_1.fq.gz and sample_name_2.fq.gz
        temp_r1_name = f"{output_prefix}_1.fq.gz"
        temp_r2_name = f"{output_prefix}_2.fq.gz"

        if os.path.exists(temp_r1_name):
            os.rename(temp_r1_name, unaligned_r1_out)
            logging.debug(f"Renamed {temp_r1_name} to {unaligned_r1_out}")
        else:
             logging.warning(f"Expected Bowtie2 output file not found: {temp_r1_name}")

        if os.path.exists(temp_r2_name):
             os.rename(temp_r2_name, unaligned_r2_out)
             logging.debug(f"Renamed {temp_r2_name} to {unaligned_r2_out}")
        else:
             logging.warning(f"Expected Bowtie2 output file not found: {temp_r2_name}")

        # --- Parse Bowtie2 log for stats ---
        total_pairs = 0
        aligned_pairs = 0
        # Updated regex to match Bowtie2's output format for total pairs
        search_pairs = re.search(r"(\d+)\s+\(.*%\) were paired", result.stderr)
        # Updated regexes to handle format like "511 (0.01%) aligned ..."
        search_aligned = re.search(r"^\s*(\d+)\s+\(.*\)\s+aligned concordantly exactly 1 time", result.stderr, re.MULTILINE)
        search_aligned_multi = re.search(r"^\s*(\d+)\s+\(.*\)\s+aligned concordantly >1 times", result.stderr, re.MULTILINE)
        # Updated regex for discordant pairs, ensuring it matches the specific indented line
        search_aligned_discord = re.search(r"^\s+(\d+)\s+\(.*\)\s+aligned discordantly 1 time", result.stderr, re.MULTILINE)

        if search_pairs:
            total_pairs = int(search_pairs.group(1))
        else:
            # Fallback: try parsing the initial read count if the 'paired' line isn't found
            # This might overestimate if there were single-end reads, but better than nothing.
            search_reads = re.search(r"(\d+) reads; of these:", result.stderr)
            if search_reads:
                total_pairs = int(search_reads.group(1))
                logging.warning(f"{sample_name}: Could not parse 'were paired' line, using total reads ({total_pairs}) as approximation for total pairs.")
            else:
                 logging.warning(f"{sample_name}: Could not parse total pairs or total reads from Bowtie2 log.")
                 total_pairs = 0 # Ensure it's zero if unparsable

        if search_aligned:
            aligned_pairs += int(search_aligned.group(1))
        if search_aligned_multi:
            aligned_pairs += int(search_aligned_multi.group(1))
        if search_aligned_discord:
            # This is an approximation, as discordantly aligned pairs might still have one mate aligned
            aligned_pairs += int(search_aligned_discord.group(1)) # Add discordant as potentially aligned

        removed_percent = 0
        if total_pairs > 0:
            # Ensure aligned_pairs doesn't exceed total_pairs due to approximations
            aligned_pairs = min(aligned_pairs, total_pairs)
            removed_percent = (aligned_pairs / total_pairs) * 100
            logging.info(f"{sample_name}: Total pairs ~={total_pairs}, Host-aligned pairs ~={aligned_pairs} ({removed_percent:.2f}% removed)")
        elif search_pairs is None: # Only log the specific warning if search_pairs failed
             # The more specific warnings about parsing are logged above
             pass

        return removed_percent, total_pairs

    except subprocess.CalledProcessError as e:
        logging.error(f"Bowtie2 failed for {sample_name}!")
        logging.error(f"Command: {' '.join(e.cmd)}")
        logging.error(f"Return code: {e.returncode}")
        logging.error(f"Stderr: {e.stderr}")
        # Write Bowtie2's stderr to a log file even on failure
        with open(log_path, 'w') as log_f:
            log_f.write(e.stderr or "No stderr captured.")
        return None, 0
    except FileNotFoundError:
        logging.error(f"Bowtie2 executable not found at '{BOWTIE2_EXE}'. Please check configuration or system PATH.")
        sys.exit(1) # Stop the script if Bowtie2 isn't found
    except Exception as e:
        logging.error(f"An unexpected error occurred during Bowtie2 run for {sample_name}: {e}")
        return None, 0


# --- Main Execution ---
if __name__ == "__main__":
    """Main function to orchestrate host removal."""
    logging.info("Starting host read removal process.")
    logging.info(f"Project base directory: {BASE_PROJECT_DIR}")
    logging.info(f"Dataset base directory: {DATASET_BASE_DIR}")
    logging.info(f"Host index: {HOST_INDEX_PATH}")
    logging.info(f"Log file: {LOG_FILE}")

    if not os.path.exists(str(HOST_INDEX_PATH) + ".1.bt2"):
         logging.error(f"Host index not found at prefix: {HOST_INDEX_PATH}")
         logging.error("Please ensure the Bowtie2 index files exist and the HOST_INDEX_PATH is set correctly (without extensions like .1.bt2).")
         sys.exit(1)

    # Determine which dataset directories to process
    if DATASET_IDS:
        datasets_to_process = [DATASET_BASE_DIR / ds_id for ds_id in DATASET_IDS]
    else:
        # Find all directories in DATASET_BASE_DIR that look like datasets (e.g., start with SG_)
        datasets_to_process = [d for d in DATASET_BASE_DIR.iterdir() if d.is_dir() and d.name.startswith("SG_")]
        if not datasets_to_process:
            logging.warning(f"No dataset directories found in {DATASET_BASE_DIR} matching 'SG_*' prefix.")

    logging.info(f"Found {len(datasets_to_process)} dataset directories to process:")
    for ds_path in datasets_to_process:
        logging.info(f"- {ds_path.name}")

    removal_stats = {} # Overall stats for final log summary
    datasets_processed_count = 0

    for dataset_dir in datasets_to_process:
        dataset_name = dataset_dir.name
        logging.info(f"--- Processing Dataset: {dataset_name} ---")

        input_fastq_dir = dataset_dir / INPUT_FASTQS_SUBDIR
        output_fastq_dir = dataset_dir / OUT_FASTQS_SUBDIR
        tsv_summary_path = dataset_dir / "about_host_read_removal.tsv" # Path for the dataset-specific TSV

        if not input_fastq_dir.is_dir():
            logging.warning(f"Input directory not found: {input_fastq_dir}. Skipping dataset {dataset_name}.")
            continue

        # Create output directory if it doesn't exist
        output_fastq_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"Output directory: {output_fastq_dir}")

        # Find FASTQ pairs
        fastq_pairs = find_fastq_pairs(str(input_fastq_dir), R1_PATTERN, R2_PATTERN)

        if not fastq_pairs:
            logging.warning(f"No valid FASTQ pairs found in {input_fastq_dir}. Skipping dataset {dataset_name}.")
            continue

        dataset_sample_stats = {} # Store stats for this dataset's TSV

        # Process each pair
        for sample_name, (r1_file, r2_file) in fastq_pairs.items():
            logging.info(f"Processing sample: {sample_name}")
            percent_removed, total_pairs = run_bowtie2(
                sample_name,
                r1_file,
                r2_file,
                HOST_INDEX_PATH,
                str(output_fastq_dir),
                NUM_THREADS,
                BOWTIE2_SENSITIVITY,
                BOWTIE2_EXE
            )
            if percent_removed is not None:
                # Log the result immediately for console visibility
                logging.info(f"Finished {sample_name}: {percent_removed:.2f}% reads removed.")
                dataset_sample_stats[sample_name] = (percent_removed, total_pairs)
                # Also add to the global stats for the final summary log
                removal_stats[f"{dataset_name}/{sample_name}"] = percent_removed
            else:
                # Log failure if needed, though run_bowtie2 already logs errors
                logging.warning(f"Processing failed for sample: {sample_name}")

        # --- Write Dataset TSV Summary ---
        if dataset_sample_stats:
            try:
                with open(tsv_summary_path, 'w') as tsv_f:
                    tsv_f.write("Sample\tPercent_Host_Reads_Removed\tTotal_Pairs\n")
                    # Sort by sample name for consistent order
                    for sample, (percent, total) in sorted(dataset_sample_stats.items()):
                        tsv_f.write(f"{sample}\t{percent:.2f}\t{total}\n")
                logging.info(f"Wrote removal summary for dataset {dataset_name} to: {tsv_summary_path}")
                datasets_processed_count += 1
            except IOError as e:
                logging.error(f"Failed to write TSV summary for {dataset_name} to {tsv_summary_path}: {e}")
        else:
            logging.warning(f"No sample statistics collected for dataset {dataset_name}, skipping TSV generation.")

    # --- Final Summary ---
    logging.info("--- Host Removal Summary ---")
    logging.info(f"Successfully processed {datasets_processed_count} dataset(s).")
    if removal_stats:
        logging.info("Overall removal percentages per sample:")
        # Sort overall stats for readability in the log
        for sample_id, percent in sorted(removal_stats.items()):
            logging.info(f"{sample_id}: {percent:.2f}% reads removed")
    else:
        logging.info("No samples were successfully processed across any dataset.")

    logging.info("Host read removal process finished.")
    logging.info(f"Overall summary log saved to: {LOG_FILE}")







