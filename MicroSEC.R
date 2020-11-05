# MicroSEC: Microhomology-induced chimeric read-originating sequence error cleaning pipeline for FFPE samples
#
# Version 1.0.6: probability calculation considering coverage depth around mutations (Nov. 3, 2020) 
# Version 1.0.5: save function modified (Nov. 2, 2020)
# Version 1.0.4: package install via github (Nov. 1, 2020)
# Version 1.0.3: parameter justification, version management with github (Oct. 29, 2020)
# Version 1.0.2: Indel length correction based on Bowtie2 penalty score (Oct. 28, 2020)
# Version 1.0.1: Correspond to the mouse genome. Homopolymer should have >=12 bases (Oct. 25, 2020)
# Version 1.0.0: Original version (Oct. 23, 2020)
# 
# Author: "Masachika Ikegami"
# 
# This pipeline is designed for filtering mutations found in formalin-fixed and paraffin-embedded (FFPE) samples.
# The MicroSEC filter utilizes a statistical analysis, and the results for mutations with less than 10 supporting reads are not reliable.
# 
# Four files are nessesary for the analysis: mutation information file, BAM file, and mutation supporting read ID information file.
#
# File 1: mutation information file
# This excel file should contain at least these contents:
#       Sample     Gene HGVS.p Mut_type Total_QV>=20   %Alt  Chr       Pos Ref Alt SimpleRepeat_TRF                     Neighborhood_sequence  Transition
# SL_1010-N6-B SLC25A24      _    1-snv          366 1.0929 chr1 108130741   C   T                N CTACCTGGAGAATGGGCCCATGTGTCCAGGTAGCAGTAAGC  C>T_t
# Total_QV>=20: The read number with total Q-value >=20. 
# SimpleRepeat_TRF: Whether the mutation locates at a simple repeat sequence or not.
# Neighborhood_sequence: [5'-20 bases] + [Alt sequence] + [3'-20 bases]
# Transition: 1-snv mutation pattern with a 3'-base. C>T_t represents CT to TT mutation. C>T_g_FFPE represents the possible FFPE artifact.
#
# File 2: BAM file
#
# File 3: mutation supporting read ID information file
# This file should contain at least these contents:
#  Chr     Pos Ref Alt                                                                                                Mut_ID     Mut
# chr1 2561609   T   A  _;ID001-1:579185f,ID004-1:1873933f;ID006-1:1131647f,ID001-1:570086f,ID008-1:1953407r,ID002-2:749570r  .;A;N#
#
# File 4: sample information tsv file  
# Seven columns are necessary.  
# [sample name] [mutation information excel file] [BAM file] [read ID information directory] [read length] [adaptor sequence] [sample type: Human or Mouse]
# PC9	./source/CCLE.xlsx	./source/Cell_line/PC9_Cell_line_Ag_TDv4.realigned.bam	./source/PC9_Cell_line	127	AGATCGGAAGAGC	Human
#
# This pipeline contains 8 filtering processes.
#
# Filter 1  : Shorter-supporting lengths distribute too short to occur (1-1 and 1-2).  
# Filter 1-1: P-values are less than the threshold_p(default: 10^(-6)).  
# Filter 1-2: The longest shorter-supporting lengths is shorter than 40% of the read length.  
# Filter 2  : Hairpin-structure induced error detection (2-1 or 2-2).  
# Filter 2-1: Palindromic sequences exist within 150 bases. 
# Filter 2-2: >=50% mutation-supporting reads contains a reverse complementary sequence of the opposite strand consisting >= 15 bases.  
# Filter 3  : 3’-/5’-supporting lengths are too densely distributed to occur (3-1, 3-2, and 3-3).  
# Filter 3-1: P-values are less than the threshold_p(default: 10^(-6)).  
# Filter 3-2: The distributions of 3’-/5’-supporting lengths are shorter than 80% of the read length.  
# Filter 3-3: <10% of bases are low quality (Quality score <18).
# Filter 4  : >=90% mutation-supporting reads are soft-clipped (after cutting adaptor sequence).  
# Filter 5  : >=20% mutations were called by chimeric reads comprising two distant regions.
# Filter 6  : Mutations locating at simple repeat sequences.
# Filter 7  : C>T_g false positive calls in FFPE samples.
# Filter 8  : Mutations locating at a >=15 homopolymer.
#
# Supporting lengths are adjusted considering small repeat sequences around the mutations.
#
# How to use
# Rscript MicroSEC.R [working/output directory] [sample information tsv file] [progress bar Y/N]
#
# Example
# Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/Sample_list.txt Y
# Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/TOP_MicroSEC_1_20/sample_info.tsv Y
# Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/TOP_MicroSEC_1_20/sample_info_test.tsv Y
#
# If you want to know the progress visually, [progress bar Y/N] should be Y.
#
# Results are saved in a excel file.
# The explatation of the results is written in detail in the second sheet of the excel file.

# load necessary packages
library(MicroSEC)

# set arguments
args = commandArgs(trailingOnly = T)
wd = args[1]
SAMPLE_LIST = args[2]
PROGRESS_BAR = args[3]

setwd(wd)
  
# initialize
MSEC = NULL
Homology_search = NULL
Mut_depth = NULL

# load sample information tsv file
SAMPLE_INFO = read.csv(SAMPLE_LIST,
                       header=FALSE,
                       stringsAsFactors=FALSE,
                       sep="\t")

for(SAMPLE in 1:dim(SAMPLE_INFO)[1]){
  SAMPLE_NAME = SAMPLE_INFO[SAMPLE,1]
  MUTATION_FILE = SAMPLE_INFO[SAMPLE,2]
  BAM_FILE = SAMPLE_INFO[SAMPLE,3]
  MUTATION_SUPPORTING_READ_LIST = SAMPLE_INFO[SAMPLE,4]
  READ_length = as.integer(SAMPLE_INFO[SAMPLE,5])
  ADAPTOR_SEQ = SAMPLE_INFO[SAMPLE,6]
  GENOME = SAMPLE_INFO[SAMPLE,7]
  
  # load mutation information
  df_mutation = fun_load_mutation(MUTATION_FILE)
  df_BAM = fun_load_BAM(BAM_FILE)
  df_mut_call = fun_load_ID(MUTATION_SUPPORTING_READ_LIST)
  
  # load genomic sequence
  genome = fun_load_genome(GENOME)
  Chr_No = fun_load_chr_no(GENOME)
  
  # analysis
  result = fun_read_check(df_mutation = df_mutation,
                          df_BAM =  df_BAM,
                          df_mut_call = df_mut_call,
                          genome = genome,
                          Chr_No = Chr_No,
                          SAMPLE_NAME = SAMPLE_NAME,
                          READ_length = READ_length,
                          ADAPTOR_SEQ = ADAPTOR_SEQ,
                          Short_Homology_search_length = 4,
                          PROGRESS_BAR = PROGRESS_BAR)
  MSEC = rbind(MSEC, result[[1]])
  Homology_search = rbind(Homology_search, result[[2]])
  Mut_depth = rbind(Mut_depth, result[[3]])
}
# search homologous sequences
MSEC = fun_homology(MSEC,
                    Homology_search,
                    Minimum_Homology_search_length = 40,
                    PROGRESS_BAR = PROGRESS_BAR)
  
# statistical analysis
MSEC = fun_summary(MSEC)
MSEC = fun_analysis(MSEC,
                    Mut_depth,
                    Short_Homology_search_length = 4,
                    Minimum_Homology_search_length = 40,
                    threshold_p = 10^(-6),
                    threshold_hairpin_ratio = 0.50,
                    threshold_soft_clip_ratio = 0.90,
                    threshold_short_length = 0.8,
                    threshold_distant_homology = 0.2,
                    threshold_low_quality_rate = 0.1,
                    Homopolymer_length = 15)
  
# save the results
fun_save(MSEC, wd)
