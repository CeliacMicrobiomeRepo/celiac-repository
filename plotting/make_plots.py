"""
Script that generates multiple plots to overview the datasets.

Plots:
1. Number of celiac samples (diagnosed or will develop) available over time
2. Celiac samples across the world (map)
3. Datasets across the world (map)
4. Plot celiac samples per sample site
5. Plot table of dataset types
6. Plot table of non-prospective sample types
7. Plot datasets by amplicon region
8. Plot datasets by sample site
9. Plot table of prospective sample types
10. Analyse and plot sex proportions
11. Analyse and plot age distributions
12. Analyse read retention over steps for 16S
13. Analyse host read removal for Shotgun
14. Analyse numbers of reads remaining after DADA2 for 16S
15. Analyse numbers of unique ASVs per sample for 16S
16. Analyse numbers of unique SGBs per sample for Shotgun
17. Analyse numbers of ASVs classified at genus and species level for 16S
18. Analyse Extraction Kit type against read counts for 16S
19. Plot samples per dataset showing significant factors
20. Plot samples per group
21. Plot Total_Num_ASVs per dataset
22. Plot Total_Num_SGBs per dataset
23. Plot age bucket proportions for celiac samples
"""


# Import libraries
import pandas as pd
import os
from plot_functions import *

# Get the directory containing this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)

# Columns: Sample_ID	Dataset_ID	SRA_Run_ID	SRA_Project_ID	Month_of_Publication	Publication_DOI	Sequencing_Type	Amplicon_Region	Num_Reads_Input	Num_Reads_Filtered	Num_Reads_DenoisedF	Num_Reads_Nonchim	Total_Pairs_Pre_Host_Removal	Percent_Host_Reads_Removed	Percent_Unclassified_Reads	Num_SGBs	Num_ASVs	Seq_Tech	DNA_Ext_Kit	DNA_Extraction_Is_Mechanical	Paired_Reads	Sample_Site	Diagnosed_Celiac	Gluten_Free_Diet	Will_Develop_Celiac	Group	Short_term_Gluten_Challenge	NCGS	Other_Autoimmune	Hookworm	Possible_Celiac	Any_Significant_Factor	Country	Age	Sex	DOID	EFO	UBERON	NCIT_Sex
# Example row: SRR5514924	16S_102_Bodkhe	SRR5514924	PRJNA385740	Feb-19	10.3389/fmicb.2019.00164	16S	V4	978413	920491	858075	266245	NA	NA	NA	NA	397	Illumina MiSeq	DNeasy Blood & Tissue Kit	False	True	duodenal	False	False	NA	HC	False	False	False	False	False	False	India	12-55	male	NA	NA	UBERON:0002114	NCIT:C20197
ALL_SAMPLES_TSV = os.path.join(PARENT_DIR, "all_samples.tsv")

# Columns: Dataset ID	Sample ID	Number of reads (nonchim)
# Example row: 16S_25_Francavilla	SRR1107516	240 reads
LOW_READ_SAMPLES_TSV = os.path.join(PARENT_DIR, "low_read_samples.tsv")

# Columns: Dataset_ID	Bioproject_ID	Record_Link	Publication_Title	Publication_Link	Month_Of_Publication	DOI	Used_In_Previous_Meta_Analysis	Lit_Search_Source	Raw_Data_Source	Essential_Metadata_Source	Sequencing_Type	Sequencing_Technology	Prospective_Study	Sample_Sites	Amplicon_Region	Forward_Primer	Reverse_Primer	DNA_Extraction_Kit	DNA_Extraction_Is_Mechanical	Read_Pairing	Trimming_Of_Reads_After_Acquisition	Bowtie2_Alignment_Sensitivity	Host_Genome_Index	MetaPhlAn_Database	Fw_Read_Trim_Position	Rv_Read_Trim_Position	ASV_Table_Length_Filter	Notes_From_Processing	Median_Num_SGBs	Median_Num_ASVs	Total_Num_ASVs	Num_ASVs_Classified_Family	Num_ASVs_Classified_Genus	Num_ASVs_Classified_Species	Age_Range	Has_Sex_Metadata	Has_Age_Metadata	Num_Samples	Num_Individuals	Num_Celiac_Samples	Num_GFD_Samples	Num_Prospective_Celiac_Samples	Longitudinal_Study	Country	Samples_With_Significant_Factors	Prospective_Studies	Shotgun_Studies	Study_Design_Description	Multiple_Publications
# Example row: 16S_136_Nobel	PRJNA778253	https://www.ncbi.nlm.nih.gov/bioproject/778253	Lack of Effect of Gluten Challenge on Fecal Microbiome in Patients With Celiac Disease and Non-Celiac Gluten Sensitivity	https://journals.lww.com/ctg/fulltext/2021/12000/lack_of_effect_of_gluten_challenge_on_fecal.5.aspx	Dec-21	10.14309/ctg.0000000000000441	FALSE	NCBI SRA	NCBI SRA	NCBI SRA	16S	Illumina MiSeq	FALSE	stool	V3-V4	TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG	GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC	MagAttract PowerSoil DNA Kit	TRUE	paired	FALSE	NA	NA	NA	0	275	437:467	-	NA	875.5	19917	18521	15499	131	34.5-76.5	FALSE	FALSE	136	25	56	55	NA	TRUE	USA	gluten challenge and NCWS	-	-	stool microbiota of HC (no GFD at any point), CD  (previously on GFD) and NCWS (previously on GFD) at 5 timepoints before, during and after a 14 day gluten challenge	-
INCLUDED_DATASETS_TSV = os.path.join(PARENT_DIR, "included_datasets.tsv")

# Columns: Publication_Title	Publication_Link	Month_Of_Publication	DOI	Record_Link	Project_ID	Lit_Search_Source	Claim_Of_Data_Availabile_Upon_Request	Data_Availability	Reason_For_Exclusion	Sequencing_Type
# Example row: High-resolution analysis of the treated coeliac disease microbiome reveals strain-level variation	https://www.tandfonline.com/doi/full/10.1080/19490976.2025.2489071	Apr-25	10.1080/19490976.2025.2489071	https://ega-archive.org/studies/EGAS50000000959	EGAS50000000959 AND EGAD00001007027	Scopus	FALSE	available on EGA	large fees to obtain celiac related metadata (see ega-archive.org)	SG
EXCLUDED_DATASETS_TSV = os.path.join(PARENT_DIR, "excluded_datasets.tsv")

# Output directory - keep plots in the plots directory
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "plots")

# Whether to skip the ASV alignments plot
SKIP_ASV_ALIGNMENTS = False

# Reference 16S gene from GTDB (E. coli 1525bp)
# >RS_GCF_030545895.1~NZ_JAUOMX010000042.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [location=567..2091] [ssu_len=1525] [contig_len=5607]
ASV_ALIGNMENT_FULL_LENGTH_16S_REFERENCE_GENE_NAME = "RS_GCF_030545895.1~NZ_JAUOMX010000042.1"
ASV_ALIGNMENT_FULL_LENGTH_16S_REFERENCE_GENE_SEQ = "ATGAAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGTCGAGGGGCAGCATGGTCTTAGCTTGCTAAGGCCGATGGCGACCGGCGCACGGGTGAGTAACACGTATCCAACCTGCCGTCTACTCTTGGACAGCCTTCTGAAAGGAAGATTAATACAAGATGGCATCATGAGTCCGCATGTTCACATGATTAAAGGTATTCCGGTAGACGATGGGGATGCGTTCCATTAGATAGTAGGCGGGGTAACGGCCCACCTAGTCTTCGATGGATAGGGGTTCTGAGAGGAAGGTCCCCCACATTGGAACTGAGACACGGTCCAAACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTCGGGTATGGATACCCGTTTGCATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGATGTTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATATCTTGAGTGCAGTTGAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCCTGCTAAGCTGCAACTGACATTGAGGCTCGAAAGTGTGGGTATCAAACAGGATTAGATACCCTGGTAGTCCACACGGTAAACGATGAATACTCGCTGTTTGCGATATACTGCAAGCGGCCAAGCGAAAGCGTTAAGTATTCCACCTGGGGAGTACGCCGGCAACGGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGAGGAACATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCCGGGCTTAAATTGCAGATGAATTACGGTGAAAGCCGTAAGCCGCAAGGCATCTGTGAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCGGCTTAAGTGCCATAACGAGCGCAACCCTTGTTGTCAGTTACTAACAGGTTCCGCTGAGGACTCTGACAAGACTGCCATCGTAAGATGTGAGGAAGGTGGGGATGACGTCAAATCAGCACGGCCCTTACGTCCGGGGCTACACACGTGTTACAATGGGGGGTACAGAGGGCCGCTACCACGCGAGTGGATGCCAATCCCCAAAACCTCTCTCAGTTCGGACTGGAGTCTGCAACCCGACTCCACGAAGCTGGATTCGCTAGTAATCGCGCATCAGCCACGGCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGGAGCCGGGGGTACCTGAAGTGCGTAACCGCGAGGAGCGCCCTAGGGTAAAACTGGTGACTGGGGCTAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGCGGCTGGAACACCTCCTTT"

# Maximum number of errors (indels + substitutions) to allow when aligning ASVs to the reference gene
# This makes the alignment robust to slight mismatches and sequencing errors
MAX_ALIGNMENT_ERRORS = 50

# If True, the script will use Will_Develop_Celiac (prospective celiac) samples in addition to Diagnosed_Celiac (non-prospective celiac) samples
USE_WILL_DEVELOP_CELIAC_SAMPLES = True



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
    # low_read_samples_df = read_tsv_safely(LOW_READ_SAMPLES_TSV)
    included_datasets_df = read_tsv_safely(INCLUDED_DATASETS_TSV)
    # excluded_datasets_df = read_tsv_safely(EXCLUDED_DATASETS_TSV)

    # If the output directory does not exist, create it
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Plot 1: Number of celiac samples over time
    celiac_samples_df = plot_celiac_samples_over_time(all_samples_df, OUTPUT_DIR, write_tsv_path="celiac_samples_over_time.tsv", include_prospective=USE_WILL_DEVELOP_CELIAC_SAMPLES)

    # Run accumulation rate analysis
    rate_results = get_celiac_accumulation_rate(celiac_samples_df, output_path=os.path.join(OUTPUT_DIR, "celiac_samples_over_time.txt"))

    # Plot 1b: Number of celiac samples over time with trend line
    plot_celiac_samples_over_time(all_samples_df, OUTPUT_DIR, plot_filename="celiac_samples_over_time_with_trend.png", trend_results=rate_results, include_prospective=USE_WILL_DEVELOP_CELIAC_SAMPLES)

    # Plot 2: Celiac samples across the world (map)
    plot_celiac_samples_across_world_map(all_samples_df, OUTPUT_DIR, write_tsv_path="celiac_samples_world_map.tsv", include_prospective=USE_WILL_DEVELOP_CELIAC_SAMPLES)

    # Plot 3: Datasets across the world (map)
    plot_datasets_across_world_map(included_datasets_df, OUTPUT_DIR, write_tsv_path="datasets_world_map.tsv")

    # Plot 4: Celiac samples per sample site
    plot_celiac_samples_per_sample_site(all_samples_df, OUTPUT_DIR, write_tsv_path="celiac_samples_per_site.tsv", include_prospective=USE_WILL_DEVELOP_CELIAC_SAMPLES)

    # Plot 5: Table of dataset types
    plot_dataset_types(included_datasets_df, OUTPUT_DIR, write_tsv_path="dataset_types.tsv")

    # Plot 6: Table of non-prospective sample types
    plot_non_prospective_sample_types(all_samples_df, OUTPUT_DIR, write_tsv_path="non_prospective_sample_types.tsv")

    # Plot 7: Datasets by amplicon region
    plot_datasets_by_amplicon_region(included_datasets_df, OUTPUT_DIR, write_tsv_path="datasets_by_amplicon_region.tsv")

    # Plot 7b: Datasets by amplicon region coloured by technology
    plot_datasets_by_amplicon_region_coloured_by_technology(included_datasets_df, OUTPUT_DIR, write_tsv_path="datasets_by_amplicon_region_coloured_by_technology.tsv")

    # Plot 8: Datasets by sample site
    plot_datasets_by_site(included_datasets_df, OUTPUT_DIR, write_tsv_path="datasets_by_site.tsv")

    # Plot 9: Table of prospective sample types
    plot_prospective_sample_types(all_samples_df, OUTPUT_DIR, write_tsv_path="prospective_sample_types.tsv")

    # Plot 10: Analyse Sex
    analyse_sex(all_samples_df, included_datasets_df, OUTPUT_DIR, write_tsv_path="sex_proportions.tsv")

    # Plot 11: Analyse Age
    analyse_age(all_samples_df, included_datasets_df, OUTPUT_DIR, write_tsv_path="age_distribution.tsv")

    # Plot 12: Analyse read retention over steps for 16S
    analyse_read_retention_16S(all_samples_df, OUTPUT_DIR, write_tsv_path="read_retention_16S.tsv")

    # Plot 13: Analyse host read removal for Shotgun
    analyse_host_read_removal_shotgun(all_samples_df, OUTPUT_DIR, write_tsv_path="host_read_removal.tsv")

    # Plot 14: Analyse numbers of reads remaining after DADA2 for 16S
    analyse_reads_remaining_after_dada2_16S(all_samples_df, OUTPUT_DIR, write_tsv_path="reads_remaining_after_dada2_16S.tsv")

    # Plot 15: Analyse numbers of unique ASVs per sample for 16S
    analyse_num_asvs_16S(all_samples_df, OUTPUT_DIR, write_tsv_path="num_asvs.tsv")

    # Plot 16: Analyse numbers of unique SGBs per sample for Shotgun
    analyse_num_sgbs_shotgun(all_samples_df, OUTPUT_DIR, write_tsv_path="num_sgbs.tsv")

    # Plot 17: Analyse numbers of ASVs classified at genus and species level for 16S
    analyse_asv_classification_16S(included_datasets_df, OUTPUT_DIR, write_tsv_path_absolute="asv_classification_16S_absolute.tsv", write_tsv_path_relative="asv_classification_16S_relative.tsv")

    # Plot 18: Analyse Extraction Kit type against read counts for 16S
    analyse_extraction_kit_type(all_samples_df, OUTPUT_DIR, write_tsv_path="reads_extraction_kit_type.tsv")

    # Plot 19: Plot samples per dataset showing significant factors
    plot_samples_per_dataset_significant_factors(all_samples_df, OUTPUT_DIR, write_tsv_path="samples_per_dataset_significant_factors.tsv")

    # Plot 20: Plot samples per group
    plot_samples_per_group(all_samples_df, OUTPUT_DIR, write_tsv_path="samples_per_group.tsv")

    # Plot 21: Plot Total_Num_ASVs per dataset
    plot_num_asvs_per_dataset(included_datasets_df, OUTPUT_DIR, write_tsv_path="num_asvs_per_dataset.tsv")

    # Plot 22: Plot Total_Num_SGBs per dataset
    plot_num_sgbs_per_dataset(included_datasets_df, OUTPUT_DIR, write_tsv_path="num_sgbs_per_dataset.tsv")

    # Plot 23: Plot age bucket proportions for celiac samples
    plot_age_bucket_proportions(all_samples_df, OUTPUT_DIR, write_tsv_path="age_bucket_proportions.tsv")
