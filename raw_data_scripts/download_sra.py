"""
Script to download raw sequencing data from NCBI SRA.

This script automates the process of downloading 16S/shotgun data,
and converting it from SRA format to FASTQ.

Usage:
1.  Configure the paths and options in the "CONFIGURATION" section below.
    - Set `DATASET_DIRS` to the list of dataset directories to process.
    - Ensure each dataset directory contains an accession list file (e.g., "SraAccList.csv").
2.  Run the script from the command line:
    `python SG_preprocessing_scripts/01_download_sra.py`
3.  The script will create subdirectories for SRA files and FASTQ files
    within each dataset directory and process the data.

Obtaining accession lists:
- These files can be downloaded from the NCBI BioProject site directly 
    - e.g. go to https://www.ncbi.nlm.nih.gov/bioproject/812940
    - click the number of SRA Experiments 
    - and press Send to > File > Format > Accession List > Create File

Example Input:
.
└── 16S_datasets
    └── 16S_102_Bodkhe
        └── SraAccList.csv

Example Output:
.
└── 16S_datasets
    └── 16S_102_Bodkhe
        ├── SraAccList.csv
        └── fastqs
            ├── SRR1039509_1.fq.gz
            ├── SRR1039509_2.fq.gz
            ├── SRR1039510_1.fq.gz
            ├── SRR1039510_2.fq.gz
            ├── SRR1039515_1.fq.gz
            ├── SRR1039515_2.fq.gz
            └── SRR1039516_1.fq.gz
            └── SRR1039516_2.fq.gz

Requirements:
-   SRA Toolkit (for `prefetch` and `fastq-dump`): 
    - https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
    - `prefetch` is used to download SRA files
    - `fastq-dump` is used to convert SRA files to FASTQ
"""


# Imports
import subprocess
import csv
from pathlib import Path
import logging
import shutil


# --- CONFIGURATION ---

# Datasets to download (must contain an SraAccList.csv file)
DATASET_DIRS = [
    # 16S datasets (on SRA) ------------
    # "16S_datasets/16S_20_Caminero/",
    # "16S_datasets/16S_20_Rawson/",
    # "16S_datasets/16S_24_Giacomin/",
    # "16S_datasets/16S_25_Francavilla/",
    # "16S_datasets/16S_27_Federica/",
    # "16S_datasets/16S_39_Olivares/",
    # "16S_datasets/16S_49_Turjeman/",
    # "16S_datasets/16S_50_Bibbo/",
    # "16S_datasets/16S_56_Laffaldano/",
    # "16S_datasets/16S_60_Shi/",
    # "16S_datasets/16S_62_Tian/",
    # "16S_datasets/16S_68_Girdhar/",
    # "16S_datasets/16S_80_Garcia/",
    # "16S_datasets/16S_96_Quagliariello/",
    # "16S_datasets/16S_99_Lionetti/",
    # "16S_datasets/16S_102_Bodkhe/",
    # "16S_datasets/16S_136_Nobel/",
    # "16S_datasets/16S_179_Verdu/",
    # "16S_datasets/16S_227_Oscarsson/",
    # "16S_datasets/16S_1211_Milletich/",

    # 16S datasets (NOT on SRA) ------------
    # "16S_datasets/16S_5_Senicar/",         # (raw data was obtained via email)
    # "16S_datasets/16S_27_Fornasaro/",      # (raw data was obtained via email)
    # "16S_datasets/16S_119_Salamon/",       # (raw data was obtained via online database: https://portalwiedzy.cm-uj.krakow.pl/info/researchdata/UJCM77a8979a493e4aacbdceefa5121abbff/)

    # Shotgun datasets (on SRA) ------------
    # "SG_datasets/SG_39_Judith/",
    # "SG_datasets/SG_80_Mouzan/"
    # "SG_datasets/SG_95_Costigan/",
    # "SG_datasets/SG_118_Leonard/",
    # "SG_datasets/SG_132_Francavilla/",
] 

# Remove intermediate files after processing
# (TEMP_ACC_LIST_FILENAME, and SRA_DIR)
REMOVE_INTERMEDIATE_FILES = True

# Renames .fastq.gz to .fq.gz (or .fastq to .fq)
RENAME_FASTQ_FILES = True

# Other constants
ACC_LIST_FILENAME = "SraAccList.csv"
TEMP_ACC_LIST_FILENAME = "acc_list.txt"
SRA_DIR = "sra_downloads"
FASTQ_DIR = "fastqs"



# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')



def setup_directories(base_dir: Path):
    """Creates necessary subdirectories within the base directory."""
    sra_path = base_dir / SRA_DIR
    fastq_path = base_dir / FASTQ_DIR
    
    sra_path.mkdir(exist_ok=True)
    fastq_path.mkdir(exist_ok=True)
    
    return sra_path, fastq_path

def read_accession_list(acc_list_file: Path):
    """Reads a list of accessions from a .txt or .csv file."""
    accessions = []
    if not acc_list_file.exists():
        raise FileNotFoundError(f"Accession list file not found: {acc_list_file}")

    if acc_list_file.suffix == ".txt":
        with open(acc_list_file, "r") as f:
            for line in f:
                accession = line.strip()
                if accession:
                    accessions.append(accession)
    elif acc_list_file.suffix == ".csv":
        with open(acc_list_file, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            try:
                next(reader)  # Skip header
            except StopIteration:
                return [] # Empty file
            for row in reader:
                if row:
                    accessions.append(row[0])
    else:
        raise ValueError(f"Unsupported accession file format: {acc_list_file.suffix}")
    
    return accessions

def download_sra_files(accessions: list, sra_dir: Path, acc_list_path: Path):
    """Downloads SRA files for a list of accessions using prefetch."""
    logging.info(f"Writing accession list to {acc_list_path}")
    with open(acc_list_path, "w") as f:
        for acc in accessions:
            f.write(f"{acc}\n")

    logging.info(f"Downloading {len(accessions)} SRA files using prefetch...")
    cmd = ["prefetch", "-O", str(sra_dir), "--option-file", str(acc_list_path)]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info("SRA files downloaded successfully.")
    except subprocess.CalledProcessError as e:
        # Check if this is just a cloud configuration warning (exit code 3)
        # but files were actually downloaded successfully
        if e.returncode == 3 and "CloudMgrGetCurrentCloud" in e.stderr:
            logging.warning("Prefetch completed with cloud configuration warnings, checking if downloads appear successful...")
            logging.warning(f"Prefetch stderr: {e.stderr}")
            
            # Verify that some files were actually downloaded
            downloaded_files = list(sra_dir.glob("*//*.sra")) + list(sra_dir.glob("*//*.sralite"))
            if downloaded_files:
                logging.info(f"Found {len(downloaded_files)} downloaded SRA files despite warnings.")
                return
            else:
                logging.error("No SRA files found - download seems to have actually failed.")
                raise
        else:
            logging.error(f"Error downloading SRA files: {e.stderr}")
            raise

def dump_sra_to_fastq(accession: str, sra_dir: Path, fastq_dir: Path):
    """Converts a single .sra file to .fastq.gz using fastq-dump."""
    sra_file_dir = sra_dir / accession
    sra_file = sra_file_dir / f"{accession}.sra"
    sra_lite_file = sra_file_dir / f"{accession}.sralite"
    
    actual_sra_path = None
    if sra_file.exists():
        actual_sra_path = sra_file
    elif sra_lite_file.exists():
        actual_sra_path = sra_lite_file
    else:
        logging.warning(f"SRA file for {accession} not found in {sra_file_dir}")
        return

    logging.info(f"Dumping {actual_sra_path} to FASTQ...")
    cmd = ["fastq-dump", "--split-3", "--gzip", "--outdir", str(fastq_dir), str(actual_sra_path)]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stderr:
            logging.warning(f"fastq-dump for {accession} produced warnings: {result.stderr}")
        logging.info(f"Successfully dumped {accession} to FASTQ.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error dumping {accession}: {e.stderr}")
        return


if __name__ == "__main__":


    # Process each dataset sequentially
    for dataset_dir_str in DATASET_DIRS:
        dataset_dir = Path(dataset_dir_str)
        logging.info(f"--- Processing dataset: {dataset_dir.name} ---")


        try:
            # 1) Setup directories
            sra_dir, fastq_dir = setup_directories(dataset_dir)
            acc_list_file = dataset_dir / ACC_LIST_FILENAME

            # 2) Read accession list
            accessions = read_accession_list(acc_list_file)
            if not accessions:
                logging.warning(f"No accessions found in {acc_list_file}. Skipping dataset.")
                continue

            # 3) Download SRA files
            temp_acc_list_path = dataset_dir / TEMP_ACC_LIST_FILENAME
            download_sra_files(accessions, sra_dir, temp_acc_list_path)
            
            # 4) Dump SRA files to FASTQ
            for i, accession in enumerate(accessions):
                logging.info(f"\n--- Processing accession {i+1}/{len(accessions)}: {accession} ---")
                dump_sra_to_fastq(accession, sra_dir, fastq_dir)

            # 5) Cleanup intermediate files
            if REMOVE_INTERMEDIATE_FILES:
                if temp_acc_list_path.exists():
                    try:
                        temp_acc_list_path.unlink()
                        logging.info(f"Removed temporary accession list: {temp_acc_list_path}")
                    except OSError as e:
                        logging.error(f"Error removing {temp_acc_list_path}: {e}")
                
                if sra_dir.exists():
                    try:
                        shutil.rmtree(sra_dir)
                        logging.info(f"Removed SRA directory: {sra_dir}")
                    except OSError as e:
                        logging.error(f"Error removing directory {sra_dir}: {e}")

            # 6) Rename FASTQ files
            # (.fastq.gz -> .fq.gz AND .fastq -> .fq)
            if RENAME_FASTQ_FILES:
                for fastq_file in list(fastq_dir.glob("*.fastq.gz")):
                    new_path = fastq_file.with_name(fastq_file.name.replace(".fastq.gz", ".fq.gz"))
                    fastq_file.rename(new_path)
                    logging.info(f"Renamed {fastq_file} to {new_path}")
                for fastq_file in list(fastq_dir.glob("*.fastq")):
                    new_path = fastq_file.with_name(fastq_file.name.replace(".fastq", ".fq"))
                    fastq_file.rename(new_path)
                    logging.info(f"Renamed {fastq_file} to {new_path}")

        except Exception as e:
            logging.error(f"An error occurred while processing {dataset_dir.name}: {e}", exc_info=True)
            continue


    logging.info("--- All datasets processed. ---") 
