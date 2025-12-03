"""
Performs taxonomic profiling of shotgun metagenomic data using MetaPhlAn 4.

This script processes paired-end FASTQ files that have undergone host read
removal. It iterates through specified dataset directories, runs MetaPhlAn on
each sample pair, and then merges the individual sample profiles into a single
taxonomic abundance table for each dataset.

Assumes:
- Input FASTQ files are paired-end, gzipped (.fq.gz), and located in a
  'host_removed_fastqs' subdirectory within each dataset directory.
- A consistent naming convention for paired reads (e.g., containing '_1' and '_2').
- MetaPhlAn (v4+), its database, Bowtie2, and the 'merge_metaphlan_tables.py'
  utility are installed and accessible in the environment/PATH.

Configuration options are set in the 'CONFIGURATION' section below.
Run this script directly from your Python IDE after running the host removal script.

Directory Structure:

Root/
├── SG_preprocessing_scripts/
│   ├── host_read_removal.py
│   ├── tax_profiling.py                 # <-- This script
│   └── convert_metaphlan_gtdb.sh
├── SG_datasets/
│   ├── <Dataset_ID_1>/                     # e.g., SG_80_Mouzan
│   │   ├── input_fastqs/
│   │   ├── host_removed_fastqs/            # Input for this script
│   │   │   ├── sample1_1.fq.gz
│   │   │   ├── sample1_2.fq.gz
│   │   │   └── ...
│   │   ├── metaphlan_outputs/              # Output directory for individual sample profiles
│   │   │   ├── sample1_profile.txt
│   │   │   ├── sample1_bowtie2.sam         # Optional intermediate file
│   │   │   └── ...
│   │   └── merged_taxonomic_profile.tsv    # Merged profile for the dataset
│   └── <Dataset_ID_2>/ ...
├── host_genome/ ...
├── metaphlan_db/                           # Directory containing the MetaPhlAn database files
│   └── mpa_vJan21_CHOCOPhlAnSGB_...        # Example DB name prefix
├── host_removal_summary.log
└── taxonomic_profiling.log                 # Log file for this script



Requirements:
- MetaPhlAn
    ```
    conda install -c bioconda metaphlan (or pip install metaphlan)
    metaphlan --version
    metaphlan --install
    ```

Note on MetaPhlAn database:
- By default, the script looks for the database in a 'metaphlan_db' directory in your project root
- Ensure the variable METAPHLAN_DB_DIR points to the correct location
- If you have multiple database versions, set METAPHLAN_DB_NAME to specify which one to use
- Database download may take a while and require several GB of disk space

Post-processing:
- After running this script, you should run the `convert_metaphlan_gtdb.sh` script
  to convert SGB-based profiles to GTDB-based profiles.
- Example command: `./convert_metaphlan_gtdb.sh ../SG_datasets/SG_39_Judith/`
- See top of script for more help.
"""



# Imports
import os
import subprocess
import glob
import logging
import sys
import pandas as pd


# ================== CONFIGURATION ==================

# --- Paths ---
# Mirroring host_read_removal.py, specify the full path to each dataset directory.
# The script will process the directories listed here.
DATASET_DIRS = [
    "SG_datasets/SG_39_Judith/",
    # "SG_datasets/SG_80_Mouzan/",
    # "SG_datasets/SG_95_Costigan/",
    # "SG_datasets/SG_118_Leonard/",
    # "SG_datasets/SG_132_Francavilla/",
]


# Subdirectory within each dataset dir containing input (host-removed) FASTQs
INPUT_FASTQS_SUBDIR = "host_removed_fastqs"
# Subdirectory within each dataset dir for MetaPhlAn outputs
METAPHLAN_OUT_SUBDIR = "metaphlan_outputs"
# Log file path
LOG_FILE = "taxonomic_profiling.log"

# --- MetaPhlAn Database ---
# Path to the DIRECTORY containing the MetaPhlAn database files (e.g., *.pkl, *.bt2)
# Example: METAPHLAN_DB_DIR = "/path/to/metaphlan_databases"
METAPHLAN_DB_DIR = "/home/haig/anaconda3/lib/python3.11/site-packages/metaphlan/metaphlan_databases" # Adjust as needed
# Database version identifier used by metaphlan (usually starts with mpa_v...)
METAPHLAN_DB_NAME = "mpa_vJan25_CHOCOPhlAnSGB_202503" # Set to specific name if needed, otherwise let MetaPhlAn choose default

# --- Tools ---
# Name of the MetaPhlAn executable (usually just 'metaphlan')
METAPHLAN_EXE = "metaphlan"
# Name of the merging utility script (usually just 'merge_metaphlan_tables.py')
MERGE_EXE = "merge_metaphlan_tables.py"
# Full path to bowtie2 executable if MetaPhlAn can't find it (usually not needed)
BOWTIE2_EXE_PATH = None # Example: "/opt/miniconda3/envs/bio/bin/bowtie2"

# --- Parameters ---
# Number of threads/processes for MetaPhlAn
NUM_THREADS = 16
# Save the intermediate bowtie2 mapping file? (Can be large)
SAVE_BOWTIE2_SAM = False
# Naming convention patterns to identify read pairs from host removal step
R1_PATTERN = "_1"
R2_PATTERN = "_2"

# ================== END CONFIGURATION ==================




# --- Logging Setup ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE, mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)

# --- Helper Functions ---
def find_fastq_pairs(fastq_dir, r1_pattern, r2_pattern):
    """Finds pairs of host-removed FASTQ files."""
    # This function is identical to the one in 02_host_read_removal.py
    # It's duplicated here for script independence, but could be refactored into a shared utility
    pairs = {}
    r1_files = glob.glob(os.path.join(fastq_dir, f"*{r1_pattern}*.fq.gz"))
    logging.info(f"Found {len(r1_files)} potential R1 files in {fastq_dir}")

    for r1_file in r1_files:
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

    logging.info(f"Identified {len(pairs)} valid FASTQ pairs for profiling.")
    return pairs

def run_metaphlan(sample_name, r1_in, r2_in, db_dir, db_name, out_dir, threads, save_sam, metaphlan_exe, bowtie2_path):
    """Runs MetaPhlAn on a single sample pair."""
    profile_output = os.path.join(out_dir, f"{sample_name}_profile.txt")
    bowtie2_output = os.path.join(out_dir, f"{sample_name}_bowtie2.sam")
    mapout_file = os.path.join(out_dir, f"{sample_name}_mapout.txt")

    # MetaPhlAn expects inputs concatenated with a comma
    input_files = f"{r1_in},{r2_in}"

    cmd = [
        metaphlan_exe,
        input_files,
        "--input_type", "fastq",
        "--nproc", str(threads),
        "-o", str(profile_output)
    ]

    # Add database path
    cmd.extend(["--db_dir", str(db_dir)])

    # Specify database name if provided, otherwise let MetaPhlAn use default within db_dir
    if db_name:
        cmd.extend(["--index", db_name]) # Note: --index parameter used in v4+

    # Add Bowtie2 path if specified
    if bowtie2_path:
        cmd.extend(["--bowtie2_exe", bowtie2_path])

    # Add mapout file, required when providing multiple (comma-separated) fastq files.
    cmd.extend(["--mapout", str(mapout_file)])

    # Always specify the bowtie2 output file path when using paired reads
    # In older versions of MetaPhlAn this was '--bowtie2out', but in v4+ it is '-s'.
    cmd.extend(["-s", str(bowtie2_output)])
    # Note: MetaPhlAn requires a SAM output path with comma-separated input,
    # even if we don't want to save the file permanently.

    logging.info(f"Running MetaPhlAn for {sample_name}...")
    logging.debug(f"Command: {' '.join(cmd)}")

    try:
        # MetaPhlAn logs progress to stderr, capture it for logging if needed
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"MetaPhlAn completed successfully for {sample_name}.")
        logging.debug(f"MetaPhlAn stdout: {result.stdout}")
        if result.stderr:
            logging.debug(f"MetaPhlAn stderr: {result.stderr}") # Contains progress

        # Delete the intermediate files if not requested
        if not save_sam:
            for fpath, ftype in [(bowtie2_output, "Bowtie2 SAM"), (mapout_file, "mapout")]:
                if os.path.exists(fpath):
                    try:
                        os.remove(fpath)
                        logging.debug(f"Removed intermediate {ftype} file: {fpath}")
                    except OSError as e:
                        logging.warning(f"Could not remove intermediate {ftype} file {fpath}: {e}")

        return profile_output # Return path to the profile file

    except subprocess.CalledProcessError as e:
        logging.error(f"MetaPhlAn failed for {sample_name}!")
        logging.error(f"Command: {' '.join(e.cmd)}")
        logging.error(f"Return code: {e.returncode}")
        logging.error(f"Stderr: {e.stderr}")
        logging.error(f"Stdout: {e.stdout}")
        return None
    except FileNotFoundError:
        logging.error(f"MetaPhlAn executable not found at '{metaphlan_exe}' or merge script not found. Please check installation and system PATH.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred during MetaPhlAn run for {sample_name}: {e}")
        return None

def merge_profiles(dataset_name, metaphlan_output_dir, merge_exe):
    """Merges individual MetaPhlAn profiles for a dataset."""
    # Save the merged file directly in the dataset directory (parent of metaphlan_output_dir)
    merged_output_file = os.path.join(os.path.dirname(metaphlan_output_dir), "merged_taxonomic_profile.tsv")
    # Use glob to find all individual profile files
    profile_files = glob.glob(os.path.join(metaphlan_output_dir, "*_profile.txt"))

    if not profile_files:
        logging.warning(f"No individual profile files (*_profile.txt) found in {metaphlan_output_dir} for dataset {dataset_name}. Skipping merge.")
        return False

    cmd = [
        merge_exe
    ]
    cmd.extend(profile_files)

    logging.info(f"Merging {len(profile_files)} MetaPhlAn profiles for dataset {dataset_name}...")
    logging.debug(f"Merge command: {merge_exe} {' '.join(profile_files)} > {merged_output_file}")

    try:
        # Redirect stdout to the output file
        with open(merged_output_file, "w") as f_out:
            result = subprocess.run(cmd, check=True, stdout=f_out, stderr=subprocess.PIPE, text=True)
        logging.info(f"Successfully merged profiles to: {merged_output_file}")
        if result.stderr:
             logging.debug(f"Merge script stderr: {result.stderr}")
        return True

    except subprocess.CalledProcessError as e:
        logging.error(f"Merging profiles failed for dataset {dataset_name}!")
        logging.error(f"Command: {' '.join(e.cmd)}")
        logging.error(f"Return code: {e.returncode}")
        logging.error(f"Stderr: {e.stderr}")
        # Clean up potentially incomplete merged file
        if os.path.exists(merged_output_file):
            os.remove(merged_output_file)
        return False
    except FileNotFoundError:
        logging.error(f"Merge script '{merge_exe}' not found. Please ensure MetaPhlAn is correctly installed.")
        return False # Don't exit the whole script, just skip merging for this dataset
    except Exception as e:
        logging.error(f"An unexpected error occurred during profile merging for {dataset_name}: {e}")
        return False

def get_unclassified_percent(profile_file):
    """
    Parses a MetaPhlAn profile file to extract the percentage of unclassified reads.
    """
    try:
        with open(profile_file, 'r') as f:
            for line in f:
                if line.startswith('UNCLASSIFIED'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        return float(parts[2])
    except IOError as e:
        logging.warning(f"Error reading profile file {profile_file}: {e}")
    return None

def generate_classification_stats(dataset_dir):
    """
    Processes a single shotgun dataset to gather MetaPhlAn unclassified percentages.
    For a given dataset it will find all sample profiles, extract the unclassified
    percentage, and write to a summary file.
    """
    metaphlan_output_dir = os.path.join(dataset_dir, 'metaphlan_outputs')
    
    if not os.path.isdir(metaphlan_output_dir):
        logging.warning(f"Skipping classification stats for {os.path.basename(dataset_dir)}: 'metaphlan_outputs' directory not found.")
        return

    profile_files = glob.glob(os.path.join(metaphlan_output_dir, '*_profile.txt'))
    
    classification_data = []
    for profile_file in profile_files:
        sample_id = os.path.basename(profile_file).replace('_profile.txt', '')
        unclassified_percent = get_unclassified_percent(profile_file)

        if unclassified_percent is not None:
            classification_data.append({
                'Sample_ID': sample_id,
                'Percent_Unclassified_Reads': unclassified_percent
            })

    if classification_data:
        df = pd.DataFrame(classification_data)
        output_path = os.path.join(dataset_dir, 'about_metaphlan_classification.tsv')
        try:
            df.to_csv(output_path, sep='\t', index=False)
            logging.info(f"Successfully wrote classification stats for {os.path.basename(dataset_dir)} to {output_path}")
        except IOError as e:
            logging.error(f"Error writing classification stats to {output_path}: {e}")
    else:
        logging.warning(f"No profile files found in {metaphlan_output_dir} for classification stats.")



if __name__ == "__main__":
    
    logging.info("Starting taxonomic profiling process with MetaPhlAn.")
    logging.info(f"Using MetaPhlAn DB directory: {METAPHLAN_DB_DIR}")
    logging.info(f"Log file: {LOG_FILE}")

    if not os.path.isdir(METAPHLAN_DB_DIR):
         logging.error(f"MetaPhlAn database directory not found: {METAPHLAN_DB_DIR}")
         logging.error("Please ensure the path is correct and the database is installed.")
         sys.exit(1)

    datasets_to_process = DATASET_DIRS
    if not datasets_to_process:
        logging.warning(f"No dataset directories specified in DATASET_DIRS. Exiting.")

    logging.info(f"Found {len(datasets_to_process)} dataset directories to process:")
    for ds_path in datasets_to_process:
        logging.info(f"- {os.path.basename(os.path.normpath(ds_path))}")

    all_datasets_successful = True

    for dataset_dir in datasets_to_process:
        dataset_name = os.path.basename(os.path.normpath(dataset_dir))
        logging.info(f"--- Processing Dataset: {dataset_name} ---")

        input_fastq_dir = os.path.join(dataset_dir, INPUT_FASTQS_SUBDIR)
        output_profile_dir = os.path.join(dataset_dir, METAPHLAN_OUT_SUBDIR)

        if not os.path.isdir(input_fastq_dir):
            logging.warning(f"Input directory not found: {input_fastq_dir}. Skipping dataset {dataset_name}.")
            continue

        # Create output directory
        os.makedirs(output_profile_dir, exist_ok=True)
        logging.info(f"Output directory: {output_profile_dir}")

        # Find host-removed FASTQ pairs
        fastq_pairs = find_fastq_pairs(input_fastq_dir, R1_PATTERN, R2_PATTERN)

        if not fastq_pairs:
            logging.warning(f"No valid host-removed FASTQ pairs found in {input_fastq_dir}. Skipping dataset {dataset_name}.")
            continue

        # Process each pair with MetaPhlAn
        profile_files_generated = []
        dataset_success = True
        for sample_name, (r1_file, r2_file) in fastq_pairs.items():
            logging.info(f"Processing sample: {sample_name}")
            profile_path = run_metaphlan(
                sample_name,
                r1_file,
                r2_file,
                METAPHLAN_DB_DIR,
                METAPHLAN_DB_NAME,
                output_profile_dir,
                NUM_THREADS,
                SAVE_BOWTIE2_SAM,
                METAPHLAN_EXE,
                BOWTIE2_EXE_PATH
            )
            if profile_path:
                profile_files_generated.append(profile_path)
            else:
                dataset_success = False # Mark dataset as failed if any sample fails
                all_datasets_successful = False

        # Merge profiles for the dataset if all samples ran successfully (or at least generated files)
        if dataset_success and profile_files_generated:
            logging.info(f"Attempting to merge profiles for dataset {dataset_name}...")
            merge_success = merge_profiles(dataset_name, output_profile_dir, MERGE_EXE)
            if not merge_success:
                all_datasets_successful = False # Mark overall process as failed if merge fails
        elif not profile_files_generated:
             logging.warning(f"No profile files were generated for dataset {dataset_name}, cannot merge.")
             all_datasets_successful = False
        else:
            logging.warning(f"Skipping profile merging for dataset {dataset_name} due to errors in sample processing.")

        # Generate classification stats for the dataset
        generate_classification_stats(dataset_dir)

    # --- Final Summary ---
    logging.info("--- Taxonomic Profiling Summary ---")
    if all_datasets_successful:
        logging.info("All specified datasets processed and profiles merged successfully.")
    else:
        logging.warning("Taxonomic profiling process completed with one or more errors. Please review the log.")

    logging.info(f"Taxonomic profiling process finished. Log saved to: {LOG_FILE}")
    

