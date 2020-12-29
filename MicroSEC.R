# MicroSEC: Microhomology-induced chimeric read-originating sequence error cleaning pipeline for FFPE samples
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
#       Sample     Gene HGVS.c HGVS.p Mut_type Total_QV>=20   %Alt  Chr       Pos Ref Alt SimpleRepeat_TRF                     Neighborhood_sequence  Transition
# SL_1010-N6-B SLC25A24   _  _    1-snv          366 1.0929 chr1 108130741   C   T                N CTACCTGGAGAATGGGCCCATGTGTCCAGGTAGCAGTAAGC  C>T_t
# Total_QV>=20: The read number with total Q-value >=20. 
# SimpleRepeat_TRF: Whether the mutation locates at a simple repeat sequence or not.
# Neighborhood_sequence: [5'-20 bases] + [Alt sequence] + [3'-20 bases]
# Sample, Mut_type, Chr, Pos, Ref, and Alt should be set exactly.  
# Gene, HGVS.c, HGVS.p, Total_QV>=20, %Alt, SimpleRepeat_TRF, and Transition can be set to any values.  
#
# File 2: BAM file
#
# File 3: mutation supporting read ID information file
# This file should contain at least these contents:
#  Chr     Pos Ref Alt                                                                                                Mut_ID     Mut
# chr1 2561609   T   A  _;ID001-1:579185f,ID004-1:1873933f;ID006-1:1131647f,ID001-1:570086f,ID008-1:1953407r,ID002-2:749570r  .;A;N#
#
# File 4: sample information tsv file  
# Seven or eight columns are necessary.  
# [sample name] [mutation information excel file] [BAM file] [read ID information directory] [read length] [adapter sequence read 1] [optional: adapter sequence read 2] [sample type: Human or Mouse]
# PC9	./source/CCLE.xlsx	./source/Cell_line/PC9_Cell_line_Ag_TDv4.realigned.bam	./source/PC9_Cell_line	127	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT Human
#
# Reference genome: Human (hg38), Mouse (mm10), hg19, hg38, or mm10
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
# Filter 4  : >=15% mutations were called by chimeric reads comprising two distant regions.
# Filter 5  : Mutations locating at simple repeat sequences.
# Filter 6  : C>T_g false positive calls in FFPE samples.
# Filter 7  : Mutations locating at a >=15 homopolymer.
# Filter 8  : >=10% low quality bases in the mutation supporting reads.
#
# Supporting lengths are adjusted considering small repeat sequences around the mutations.
#
# How to use
# Style 1
# Rscript MicroSEC.R [working/output directory] [sample information tsv file] [progress bar Y/N]
# 
# Style 2 (for pipeline use)
# Rscript MicroSEC.R [output filename] [sample name] [mutation information tsv.gz] [BAM file] [read ID information directory] [read length] [adapter sequence read 1] [adapter sequence read 2] [sample type: Human or Mouse]
#
# Example
# Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/Sample_list_test.txt N
# Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/TOP_MicroSEC_1_20/sample_info_test.tsv N
# Rscript MicroSEC.R /mnt/result/post_filter/SAMPLE.gz /mnt/result/mutation/SAMPLE.gz /mnt/result/mutation/SAMPLE.gz /mnt/result/BAM/SAMPLE.bam /mnt/result/ID 150 	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT hg38
#
# If you want to know the progress visually, [progress bar Y/N] should be Y.
#
# Results are saved in a excel file.
# The explatation of the results is written in detail in the second sheet of the excel file.

# load necessary packages
library(MicroSEC)

# set arguments
args = commandArgs(trailingOnly = T)

if (args[3] == "N" | args[3] == "Y") {
  wd = args[1]
  sample_list = args[2]
  progress_bar = args[3]
  
  setwd(wd)
  
  # load sample information tsv file
  sample_info = read.csv(sample_list,
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         sep = "\t")
  
  # initialize
  msec = NULL
  homology_search = NULL
  mut_depth = NULL
  
  for (sample in seq_len(dim(sample_info)[1])) {
    sample_name = sample_info[sample, 1]
    mutation_file = sample_info[sample, 2]
    bam_file = sample_info[sample, 3]
    read_list = sample_info[sample, 4]
    read_length = as.integer(sample_info[sample, 5])
    adapter_1 = sample_info[sample, 6]
    if (sample_info[sample, 7] %in%
        c("Human", "Mouse", "hg19", "hg38", "mm10")) {
      adapter_2 = adapter_1
      organism = sample_info[sample, 7]
    } else{
      adapter_2 = sample_info[sample, 7]
      organism = sample_info[sample, 8]
    }
    
    # load mutation information
    df_mutation = fun_load_mutation(mutation_file, sample_name)
    df_bam = fun_load_bam(bam_file)
    df_mut_call = fun_load_id(read_list)
    
    # load genomic sequence
    ref_genome = fun_load_genome(organism)
    chr_no = fun_load_chr_no(organism)
    
    # analysis
    result = fun_read_check(df_mutation = df_mutation,
                            df_bam =  df_bam,
                            df_mut_call = df_mut_call,
                            ref_genome = ref_genome,
                            sample_name = sample_name,
                            read_length = read_length,
                            adapter_1 = adapter_1,
                            adapter_2 = adapter_2,
                            short_homology_search_length = 4,
                            progress_bar = progress_bar)
    msec = rbind(msec, result[[1]])
    homology_search = rbind(homology_search, result[[2]])
    mut_depth = rbind(mut_depth, result[[3]])
  }
  # search homologous sequences
  msec = fun_homology(msec,
                      homology_search,
                      min_homology_search = 40,
                      ref_genome,
                      chr_no,
                      progress_bar = progress_bar)
  
  # statistical analysis
  msec = fun_summary(msec)
  msec = fun_analysis(msec,
                      mut_depth,
                      short_homology_search_length = 4,
                      min_homology_search = 40,
                      threshold_p = 10 ^ (-6),
                      threshold_hairpin_ratio = 0.50,
                      threshold_short_length = 0.75,
                      threshold_distant_homology = 0.15,
                      threshold_low_quality_rate = 0.1,
                      homopolymer_length = 15)
  
  # save the results
  fun_save(msec, sample_info[1,1], wd)
}else {
  output = args[1]
  sample_name = args[2]
  mutation_file = args[3]
  bam_file = args[4]
  read_list = args[5]
  read_length = as.integer(args[6])
  adapter_1 = args[7]
  adapter_2 = args[8]
  organism = args[9]
  progress_bar = "N"
  
  # load mutation information
  df_mutation = fun_load_mutation_gz(mutation_file)
  df_bam = fun_load_bam(bam_file)
  df_mut_call = fun_load_id(read_list)
  
  # load genomic sequence
  ref_genome = fun_load_genome(organism)
  chr_no = fun_load_chr_no(organism)
  
  # analysis
  result = fun_read_check(df_mutation = df_mutation,
                          df_bam =  df_bam,
                          df_mut_call = df_mut_call,
                          ref_genome = ref_genome,
                          sample_name = sample_name,
                          read_length = read_length,
                          adapter_1 = adapter_1,
                          adapter_2 = adapter_2,
                          short_homology_search_length = 4,
                          progress_bar = progress_bar)
  msec = result[[1]]
  homology_search = result[[2]]
  mut_depth = result[[3]]
  
  # search homologous sequences
  msec = fun_homology(msec,
                      homology_search,
                      min_homology_search = 40,
                      ref_genome,
                      chr_no,
                      progress_bar = progress_bar)
  
  # statistical analysis
  msec = fun_summary(msec)
  msec = fun_analysis(msec,
                      mut_depth,
                      short_homology_search_length = 4,
                      min_homology_search = 40,
                      threshold_p = 10 ^ (-6),
                      threshold_hairpin_ratio = 0.50,
                      threshold_short_length = 0.75,
                      threshold_distant_homology = 0.15,
                      threshold_low_quality_rate = 0.1,
                      homopolymer_length = 15)
  
  # save the results
  fun_save_gz(msec, output)
}

