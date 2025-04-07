"""
Script that generates multiple plots to overview the datasets.

Plots:
1. Number of diagnosed celiac samples available over time
2. Celiac samples across the world (map)
3. Plot celiac samples per sample site
4. Plot table of dataset types
5. Plot table of sample types
6. Plot datasets by amplicon region
"""


# Import libraries
import pandas as pd
import os
from plot_functions import *

# Get the directory containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)

# Columns: Sample_ID	Dataset_ID	SRA_Run_ID	SRA_Project_ID	Month_of_Publication	Publication_DOI	Sequencing_Type	Amplicon_Region	Seq_Depth_Filtered	Seq_Tech	DNA_Ext_Kit	Paired_Reads	Sample_Site	Diagnosed_Celiac	Gluten_Free_Diet	Group	Will_Develop_Celiac	Group_Prospective_Study	Short_term_Gluten_Challenge	NCGS	Other_Autoimmune	Hookworm	Possible_Celiac	Any_Significant_Factor	Country	Age	Sex							
# Example row: SRR5514924	16S_102_Bodkhe	SRR5514924	PRJNA385740	Feb-19	10.3389/fmicb.2019.00164	16S	V4	920491	Illumina MiSeq	DNeasy Blood & Tissue Kit	TRUE	duodenal	FALSE	FALSE	HC	NA	NA	FALSE	FALSE	FALSE	FALSE	FALSE	FALSE	India	12-55	male							
ALL_SAMPLES_TSV = os.path.join(PARENT_DIR, "all_samples.tsv")

# Columns: Dataset ID	Sample ID	Number of reads (nonchim)
# Example row: 16S_25_Francavilla	SRR1107516	240 reads
LOW_READ_SAMPLES_TSV = os.path.join(PARENT_DIR, "low_read_samples.tsv")

# Columns: dataset_ID	bioproject_ID	record link	publication title	publication link	month of publication	DOI	used in previous meta-analysis	Lit. search source	Data source	sequencing type	sequencing technology	prospective study	sample site(s)	amplicon region	V1	V2	V3	V4	V5	V6	forward primer	reverse primer	DNA extraction kit	read pairing	trimming of reads (after acquisition)	fw read trim position	rv read trim position	ASV table length filter	notes from processing	age range	num samples (processed and with metadata)	num individuals (processed and with metadata)	num celiac samples (processed and with metadata)	num GFD samples (processed and with metadata)	num prospective celiac samples (processed and with metadata)	longitudinal study	country	samples with significant factors	V4 stool	stool class imbalance	V4 duodenal	duodenal class imbalance	core V4 dataset	non-V4 studies	non-stool non-duodenal studies	prospective studies	shotgun studies	study design description	
# Example row: 16S_60_Shi	PRJNA890948	https://www.ncbi.nlm.nih.gov/bioproject/890948	Characteristics of gut microbiota and fecal metabolomes in patients with celiac disease in Northwest China	https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.1020977/full	Nov-22	10.3389/fmicb.2022.1020977	FALSE	NCBI SRA	NCBI SRA	16S	Illumina NovaSeq 6000	FALSE	stool	V3-V4			TRUE	TRUE			341F (5?-CCTACGGGNGGCWGCAG-3?)	805R (5?-GACTACHVGGGTATCTAATCC-3?)	FastDNA Spin Kit for Soil	paired	FALSE	0	0	433:467		19-65	60	60	30	0	NA	FALSE	China	-	TRUE	-	-	-	TRUE	-	-	-	-	stool microbiota of ACD and HC	
INCLUDED_DATASETS_TSV = os.path.join(PARENT_DIR, "included_datasets.tsv")

# Columns: Publication title	Publication link	Month of publication	DOI	record link	project_ID	Lit. search source	Data availability	Reason for exclusion	Sequencing type																					
# Example row: Effects of a low FODMAP diet on gut microbiota in individuals with treated coeliac disease having persistent gastrointestinal symptoms   a randomised controlled trial	https://www.cambridge.org/core/journals/british-journal-of-nutrition/article/effects-of-a-low-fodmap-diet-on-gut-microbiota-in-individuals-with-treated-coeliac-disease-having-persistent-gastrointestinal-symptoms-a-randomised-controlled-trial/BC777B9C41E45392E54CB8148BB6E3E9	Jun-23	10.1017/S0007114523001253	NA	NA	Scopus	unavailable	no email response from corresponding author(s)	16S																					
EXCLUDED_DATASETS_TSV = os.path.join(PARENT_DIR, "excluded_datasets.tsv")

# Output directory - keep plots in the plots directory
OUTPUT_DIR = SCRIPT_DIR


def read_tsv_safely(filepath):
    """Read TSV file with error handling for encoding issues."""
    try:
        # Try UTF-8 first
        return pd.read_csv(filepath, sep="\t")
    except UnicodeDecodeError:
        # If UTF-8 fails, try with 'latin-1' encoding which can handle most special characters
        return pd.read_csv(filepath, sep="\t", encoding='latin-1')


if __name__ == "__main__":

    # Load all TSV files
    all_samples_df = read_tsv_safely(ALL_SAMPLES_TSV)
    low_read_samples_df = read_tsv_safely(LOW_READ_SAMPLES_TSV)
    included_datasets_df = read_tsv_safely(INCLUDED_DATASETS_TSV)
    excluded_datasets_df = read_tsv_safely(EXCLUDED_DATASETS_TSV)

    # If the output directory does not exist, create it
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Plot 1: Number of celiac samples over time
    plot_celiac_samples_over_time(all_samples_df, OUTPUT_DIR)

    # Plot 2: Celiac samples across the world (map)
    plot_celiac_samples_across_world_map(all_samples_df, OUTPUT_DIR)

    # Plot 3: Celiac samples per sample site
    plot_celiac_samples_per_sample_site(all_samples_df, OUTPUT_DIR)

    # Plot 4: Table of dataset types
    plot_dataset_types(included_datasets_df, OUTPUT_DIR)

    # Plot 5: Table of sample types
    plot_sample_types(all_samples_df, OUTPUT_DIR)

    # Plot 6: Datasets by amplicon region
    plot_datasets_by_amplicon_region(included_datasets_df, OUTPUT_DIR)





