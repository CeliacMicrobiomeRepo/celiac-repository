"""
Script to trim raw sequencing data.

This script automates the process of trimming 16S/shotgun FASTQ
data using Trimmomatic and/or Cutadapt.

Usage:
1.  Obtain the raw FASTQ files from SRA or other source.
    - If they are on SRA, run `01_download_sra.py` to download the data.
    - Else, manually place them in the `fastqs` directory.
    - See Example Input below.
2.  Configure the paths and options in the "CONFIGURATION" section below.
    - Set `DATASET_DIRS` to the list of dataset directories to process.
    - Ensure the `fastqs` directory exists within each dataset directory and
      contains the raw FASTQ files.
    - Set other options to do with file name and trimming (see below).
3.  Run the script from the command line:
    `python SG_preprocessing_scripts/02_trim_reads.py`
4.  The script will create a subdirectory for trimmed FASTQ files. If the
    directory already contains files, you will be prompted to overwrite it.


Example Input:
.
└── 16S_datasets
    └── 16S_102_Bodkhe
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
- Trimmomatic: 
    - http://www.usadellab.org/cms/?page=trimmomatic
    - Download the Trimmomatic jar file and place it in the `Trimmomatic-0.39` directory.
- Cutadapt: 
    - https://cutadapt.readthedocs.io/en/stable/installation.html
    - Install Cutadapt for CLI
"""


# Imports
import subprocess
from pathlib import Path
import logging
import shutil


# --- CONFIGURATION ---

# General Options
DATASET_DIRS = [
    # 16S datasets (on SRA) ------------
    # "16S_datasets/16S_20_Caminero/",       # Paired-end trimmed with Cutadapt (see below)
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
    # "16S_datasets/16S_102_Bodkhe/",       # Paired-end trimmed with Trimmomatic (TruSeq3-PE-2.fa)
    # "16S_datasets/16S_136_Nobel/",
    # "16S_datasets/16S_179_Verdu/",        # Paired-end trimmed with Cutadapt (see below)
    # "16S_datasets/16S_227_Oscarsson/",
    # "16S_datasets/16S_1211_Milletich/",

    # 16S datasets (NOT on SRA) ------------
    # These would need their fastq files placed manually in the `fastqs` directory
    # "16S_datasets/16S_5_Senicar/",
    # "16S_datasets/16S_27_Fornasaro/",
    # "16S_datasets/16S_119_Salamon/",

    # Shotgun datasets (on SRA) ------------
    # "SG_datasets/SG_39_Judith/",         # ???
    # "SG_datasets/SG_80_Mouzan/"          # Paired-end trimmed with Trimmomatic (NexteraPE-PE.fa)
    # "SG_datasets/SG_95_Costigan/",       # Paired-end trimmed with Trimmomatic (NexteraPE-PE.fa)
    # "SG_datasets/SG_118_Leonard/",       # ???
    # "SG_datasets/SG_132_Francavilla/",
]


# File Name Options ------
# Input directory
FASTQ_DIR = "fastqs"
# Input files
FW_READ_SUFFIX = "_1.fq.gz"
RV_READ_SUFFIX = "_2.fq.gz"
# Output directory
TRIMMED_DIR = "trimmed"


# Trimmomatic Options ------
TRIM_WITH_TRIMMOMATIC = False
# Path to Trimmomatic jar file
TRIMMOMATIC_PATH = "Trimmomatic-0.39/trimmomatic-0.39.jar"
# Select the trimmomatic adapter file
# "TruSeq2 is used in GAII machines and TruSeq3 is used by HiSeq and MiSeq machines"
TRIMMOMATIC_ADAPTERS = "Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
# Number of threads to use
THREADS = 16
# Delete unpaired reads?
DELETE_UNPAIRED = True


# Cutadapt Options ------
TRIM_WITH_CUTADAPT = True
# Used for 16S_20_Caminero and 16S_179_Verdu
FWD_CUTADAPT_ADAPTERS = [
    ("GTATTACCGCGGCTGCTGG", 'front')
]
REV_CUTADAPT_ADAPTERS = [
    ("CCTACGGGAGGCAGCAG", 'front')
]


# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def analyse_fastq_directory(fastq_dir: Path, fwd_suffix: str, rev_suffix: str):
    """
    Analyses a directory of FASTQ files to determine samples and read type.

    This function identifies samples based on forward read suffixes, checks for
    corresponding reverse reads, and ensures consistency (all paired-end or all
    single-end) across the dataset.

    Args:
        fastq_dir: Path to the directory containing FASTQ files.
        fwd_suffix: The filename suffix for forward reads (e.g., "_1.fq.gz").
        rev_suffix: The filename suffix for reverse reads (e.g., "_2.fq.gz").

    Returns:
        A tuple containing:
        - A list of dictionaries, where each dictionary represents a sample
          and contains its 'accession', 'fwd_path', and 'rev_path' (or None).
        - A boolean indicating if the dataset is paired-end.

    Raises:
        FileNotFoundError: If the fastq_dir does not exist.
        ValueError: If a mix of single-end and paired-end reads is found.
    """
    if not fastq_dir.is_dir():
        raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")

    # Discover samples based on forward read suffix
    fwd_files = sorted(list(fastq_dir.glob(f"*{fwd_suffix}")))
    if not fwd_files:
        return [], False

    samples = []
    is_paired_list = []
    processed_files = set()

    for fwd_path in fwd_files:
        accession = fwd_path.name.replace(fwd_suffix, "")
        rev_path = fastq_dir / f"{accession}{rev_suffix}"

        sample_info = {'accession': accession, 'fwd_path': fwd_path}
        processed_files.add(fwd_path)

        if rev_path.exists():
            sample_info['rev_path'] = rev_path
            is_paired_list.append(True)
            processed_files.add(rev_path)
        else:
            sample_info['rev_path'] = None
            is_paired_list.append(False)
        
        samples.append(sample_info)

    # Verify that all samples are consistently paired or single
    if any(is_paired_list) and not all(is_paired_list):
        raise ValueError("Inconsistent read types: mixture of single-end and paired-end samples found.")

    # Warn about any FASTQ files that don't match the expected naming scheme
    all_fastq_files = set(fastq_dir.glob("*.f*q.gz"))
    unmatched_files = all_fastq_files - processed_files
    if unmatched_files:
        logging.warning(f"Found files that do not match expected suffixes and will be ignored: "
                        f"{[str(f.name) for f in unmatched_files]}")

    is_paired = all(is_paired_list) if is_paired_list else False
    return samples, is_paired


def trim_with_trimmomatic(accession: str, fwd_in: Path, rev_in: Path, out_dir: Path, is_paired: bool):
    """Trims reads using Trimmomatic."""
    logging.info(f"Trimming {accession} with Trimmomatic...")
    
    illumina_clip = f"ILLUMINACLIP:{TRIMMOMATIC_ADAPTERS}:2:30:10:4:true"
    
    if is_paired:
        fwd_out_paired = out_dir / f"{accession}{FW_READ_SUFFIX}"
        fwd_out_unpaired = out_dir / f"{accession}_1_unpaired.fq.gz"
        rev_out_paired = out_dir / f"{accession}{RV_READ_SUFFIX}"
        rev_out_unpaired = out_dir / f"{accession}_2_unpaired.fq.gz"

        cmd = [
            "java", "-jar", TRIMMOMATIC_PATH, "PE",
            "-threads", str(THREADS),
            str(fwd_in), str(rev_in),
            str(fwd_out_paired), str(fwd_out_unpaired),
            str(rev_out_paired), str(rev_out_unpaired),
            illumina_clip
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.info(f"Trimmomatic finished for {accession}.")

            if DELETE_UNPAIRED:
                logging.info(f"Deleting unpaired read files for {accession}.")
                fwd_out_unpaired.unlink(missing_ok=True)
                rev_out_unpaired.unlink(missing_ok=True)

            return fwd_out_paired, rev_out_paired
        except subprocess.CalledProcessError as e:
            logging.error(f"Trimmomatic failed for {accession}: {e.stderr}")
            return None, None

    else:  # Single-end
        fwd_out = out_dir / f"{accession}{FW_READ_SUFFIX}"
        cmd = [
            "java", "-jar", TRIMMOMATIC_PATH, "SE",
            "-threads", str(THREADS),
            str(fwd_in),
            str(fwd_out),
            illumina_clip
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.info(f"Trimmomatic finished for {accession}.")
            return fwd_out, None
        except subprocess.CalledProcessError as e:
            logging.error(f"Trimmomatic failed for {accession}: {e.stderr}")
            return None, None


def trim_with_cutadapt(accession: str, fwd_in: Path, rev_in: Path, out_dir: Path, is_paired: bool):
    """Trims reads using Cutadapt."""
    logging.info(f"Trimming {accession} with Cutadapt...")
    
    fwd_out = out_dir / f"{accession}{FW_READ_SUFFIX}"
    rev_out = out_dir / f"{accession}{RV_READ_SUFFIX}" if is_paired else None

    cmd = ["cutadapt", "--minimum-length=10", "-n=10"]
    
    if is_paired:
        for adapter, pos in FWD_CUTADAPT_ADAPTERS:
            if pos == 'front':
                cmd.extend(['-g', f'{adapter};e=0.25'])
        for adapter, pos in REV_CUTADAPT_ADAPTERS:
            if pos == 'front':
                cmd.extend(['-G', f'{adapter};e=0.25'])
        cmd.extend(['-o', str(fwd_out), '-p', str(rev_out), str(fwd_in), str(rev_in)])
    else:  # Single-end
        for adapter, pos in FWD_CUTADAPT_ADAPTERS:
            if pos == 'front':
                cmd.extend(['-g', f'{adapter};e=0.25'])
        cmd.extend(['-o', str(fwd_out), str(fwd_in)])

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"Cutadapt finished for {accession}.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Cutadapt failed for {accession}: {e.stderr}")
        return None, None

    return fwd_out, rev_out


if __name__ == "__main__":

    # --- 1. Validate Configuration ---
    if TRIM_WITH_TRIMMOMATIC and TRIM_WITH_CUTADAPT:
        logging.error("Only one of Trimmomatic or Cutadapt can be used at a time. Exiting.")
        exit()
    if not (TRIM_WITH_TRIMMOMATIC or TRIM_WITH_CUTADAPT):
        logging.error("No trimming method selected. Set TRIM_WITH_TRIMMOMATIC or TRIM_WITH_CUTADAPT to True.")
        exit()


    # --- 2. Process Each Dataset ---
    for dataset_dir_str in DATASET_DIRS:
        dataset_dir = Path(dataset_dir_str)
        logging.info(f"--- Processing dataset: {dataset_dir.name} ---")

        try:
            # --- 2a. Analyse FASTQ directory and discover samples ---
            fastq_dir = dataset_dir / FASTQ_DIR
            samples, is_paired = analyse_fastq_directory(fastq_dir, FW_READ_SUFFIX, RV_READ_SUFFIX)

            if not samples:
                logging.warning(f"No valid FASTQ files found in {fastq_dir} matching the suffixes. Skipping dataset.")
                continue
            read_type = "paired-end" if is_paired else "single-end"
            logging.info(f"Found {len(samples)} samples, detected as {read_type}.")

            # --- 2b. Set up output directory, prompting for overwrite if needed ---
            trimmed_dir = dataset_dir / TRIMMED_DIR
            if trimmed_dir.exists() and any(trimmed_dir.iterdir()):
                response = input(f"Output directory '{trimmed_dir}' is not empty. "
                                 f"Press ENTER to remove its contents and continue, or 's' to skip: ")
                if response.lower() == 's':
                    logging.warning(f"Skipping dataset {dataset_dir.name} by user request.")
                    continue
                logging.info(f"Clearing output directory: {trimmed_dir}")
                shutil.rmtree(trimmed_dir)
            
            trimmed_dir.mkdir(exist_ok=True)

            # --- 2c. Iterate through samples and trim ---
            for i, sample in enumerate(samples):
                accession = sample['accession']
                logging.info(f"--- Processing sample {i+1}/{len(samples)}: {accession} ---")

                current_fwd = sample['fwd_path']
                current_rev = sample['rev_path']
                
                if not current_fwd.exists():
                    logging.error(f"Forward read for {accession} not found at {current_fwd}. Skipping.")
                    continue
                if is_paired and not current_rev.exists():
                    logging.error(f"Reverse read for {accession} not found at {current_rev}, but in paired-end mode. Skipping.")
                    continue
                
                if TRIM_WITH_TRIMMOMATIC:
                    fwd_after_trim, rev_after_trim = trim_with_trimmomatic(accession, current_fwd, current_rev, trimmed_dir, is_paired)
                    if not fwd_after_trim:
                        logging.warning(f"Trimmomatic failed for {accession}, skipping.")
                        continue

                elif TRIM_WITH_CUTADAPT:
                    fwd_after_cutadapt, rev_after_cutadapt = trim_with_cutadapt(accession, current_fwd, current_rev, trimmed_dir, is_paired)
                    if not fwd_after_cutadapt:
                        logging.warning(f"Cutadapt failed for {accession}, skipping.")
                        continue
        
        except (FileNotFoundError, ValueError) as e:
            logging.error(f"A validation error occurred while processing {dataset_dir.name}: {e}", exc_info=True)
            continue
        except Exception as e:
            logging.error(f"An unexpected error occurred while processing {dataset_dir.name}: {e}", exc_info=True)
            continue

    logging.info("--- All datasets processed. ---") 

