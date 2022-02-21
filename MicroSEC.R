# MicroSEC: Microhomology-induced chimeric read-originating sequence error cleaning pipeline for FFPE samples
#
# Author: "Masachika Ikegami"
# 
# This pipeline is designed for filtering mutations found in formalin-fixed and paraffin-embedded (FFPE) samples.
# The MicroSEC filter utilizes a statistical analysis, and the results for mutations with less than 10 supporting reads are not reliable.
# 
# Two files are necessary for the analysis: mutation information file, BAM file
# An additional file is preferable: mutation supporting read ID information file.  
# A sample information tsv file is mandatory if multiple samples are processed simultaneously.  
#
# File 1: mutation information file (mandatory)
# This excel file should contain at least these columns:
#       Sample     Gene HGVS.c HGVS.p Mut_type Total_QV>=20   %Alt  Chr       Pos Ref Alt SimpleRepeat_TRF                     Neighborhood_sequence  Transition
# SL_1010-N6-B SLC25A24   _  _    1-snv          366 1.0929 chr1 108130741   C   T                N CTACCTGGAGAATGGGCCCATGTGTCCAGGTAGCAGTAAGC  C>T_t
# Total_QV>=20: The read number with total Q-value >=20. 
# SimpleRepeat_TRF: Whether the mutation locates at a simple repeat sequence or not.
# Neighborhood_sequence: [5'-20 bases] + [Alt sequence] + [3'-20 bases]
# Sample, Mut_type, Chr, Pos, Ref, and Alt should be set exactly.  
# Gene, HGVS.c, HGVS.p, Total_QV>=20, %Alt, SimpleRepeat_TRF, and Transition can be set to any values.  
# If you do not know the Neighborhood_sequence, enter "-".
#
# File 2: BAM file (mandatory)
#
# File 3: mutation supporting read ID information file (optional)
# This file should contain at least these contents:
#  Chr     Pos Ref Alt                                                                                                Mut_ID     Mut
# chr1 2561609   T   A  _;ID001-1:579185f,ID004-1:1873933f;ID006-1:1131647f,ID001-1:570086f,ID008-1:1953407r,ID002-2:749570r  .;A;N#
#
# File 4: sample information tsv file  (mandatory, if multiple samples are processed in a batch)
# Six to eight columns are necessary (without column names).
# Optional columns can be deleted if they are not applicable.
# [sample name] [mutation information excel file] [BAM file] [optional: read ID information directory] [read length] [adapter sequence read 1] [optional: adapter sequence read 2] [sample type: Human or Mouse]
# PC9	./source/CCLE.xlsx	./source/Cell_line/PC9_Cell_line_Ag_TDv4.realigned.bam	./source/PC9_Cell_line	127	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT Human
# A375	./source/CCLE.xlsx	./source/Cell_line/A375_Cell_line_Ag_TDv4.realigned.bam	127	 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT Human

# Reference genome: Human (hg38 or hg19) or Mouse (mm10)
#
# This pipeline contains 8 filtering processes.
#
# Filter 1  : Shorter-supporting lengths distribute too unevenly to occur (1-1 and 1-2).  
# Filter 1-1: P-values are less than the threshold_p(default: 10^(-6)).  
# Filter 1-2: The shorter-supporting lengths distributed over less than 75% of the read length.  
# Filter 2  : Hairpin-structure induced error detection (2-1 and 2-2).  
# Filter 2-1: Palindromic sequences exist within 150 bases. 
# Filter 2-2: >=50% mutation-supporting reads contains a reverse complementary sequence of the opposite strand consisting >= 15 bases.  
# Filter 3  : 3’-/5’-supporting lengths are too unevenly distributed to occur (3-1 and 3-3).  
# Filter 3-1: P-values are less than the threshold_p(default: 10^(-6)).  
# Filter 3-2: The distributions of 3’-/5’-supporting lengths are within 75% of the read length.  
# Filter 4  : >=15% mutations were called by chimeric reads comprising two distant regions.
# Filter 5  : Mutations locating at simple repeat sequences.
# Filter 6  : C>T_g false positive calls in FFPE samples.
# Filter 7  : Mutations locating at a >=15 homopolymer.
# Filter 8  : >=10% low quality bases (Quality score <18) in the mutation supporting reads.
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
# Rscript MicroSEC.R /mnt/result/post_filter/SAMPLE.gz SAMPLE /mnt/result/mutation/SAMPLE.tsv.gz /mnt/result/BAM/SAMPLE.bam /mnt/result/ID 150 	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT hg38
#
# If you want to know the progress visually, [progress bar Y/N] should be Y.
#
# Results are saved in a excel file.
# The explatation of the results is written in detail in the second sheet of the excel file.

# load necessary packages
library(MicroSEC)

# set arguments
args <- commandArgs(trailingOnly = T)

if (args[3] == "N" | args[3] == "Y") {
  wd <- args[1]
  sample_list <- args[2]
  progress_bar <- args[3]
  
  setwd(wd)
  
  # load sample information tsv file
  sample_info <- read.csv(sample_list,
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         sep = "\t")
  
  # initialize
  msec <- NULL
  homology_search <- NULL
  mut_depth <- NULL
  
  for (sample in seq_len(dim(sample_info)[1])) {
    sample_name = sample_info[sample, 1]
    mutation_file = sample_info[sample, 2]
    bam_file = sample_info[sample, 3]
    if (is.character(sample_info[sample, 4]) == TRUE) {
      list_exist <- TRUE
      read_list <- sample_info[sample, 4]
      read_length <- as.integer(sample_info[sample, 5])
      adapter_1 <- sample_info[sample, 6]
      if (sample_info[sample, 7] %in%
          c("Human", "Mouse", "hg19", "hg38", "mm10")) {
        adapter_2 <- adapter_1
        organism <- sample_info[sample, 7]
      } else{
        adapter_2 <- sample_info[sample, 7]
        organism <- sample_info[sample, 8]
      }
    } else {
      list_exist <- FALSE
      read_length <- as.integer(sample_info[sample, 4])
      adapter_1 <- sample_info[sample, 5]
      if (sample_info[sample, 6] %in%
          c("Human", "Mouse", "hg19", "hg38", "mm10")) {
        adapter_2 <- adapter_1
        organism <- sample_info[sample, 6]
      } else{
        adapter_2 <- sample_info[sample, 6]
        organism <- sample_info[sample, 7]
      }
    }

    bam_file_bai = paste(bam_file, ".bai", sep="")
    bam_file_slim = paste(bam_file, ".SLIM", sep="")
    if (!file.exists(bam_file_bai) & !file.exists(bam_file_slim)) {
      print("Sorting a BAM file...")
      bam_file_sort = paste(bam_file, "_sort.bam", sep="")
      syscom = paste("samtools sort -@ 4 -o ",
                     bam_file_sort,
                     " ",
                     bam_file,
                     sep="")
      system(syscom)
      syscom = paste("samtools index ",
                     bam_file_sort,
                     sep="")
      system(syscom)
      bam_file = bam_file_sort
    }
    bam_file_tmp1 = paste(bam_file, ".tmp1", sep="")
    bam_file_tmp2 = paste(bam_file, ".tmp2", sep="")
    bam_file_slim = paste(bam_file, ".SLIM", sep="")

    # load genomic sequence
    fun_load_genome(organism) # ref_genome
    fun_load_chr_no(organism) # chr_no
    
    # load mutation information
    fun_load_mutation(mutation_file, sample_name) # df_mutation
    df_mutation = df_mutation[order(df_mutation$Chr, df_mutation$Pos),]
    sep_new = TRUE
    continuous = FALSE
    chr_last = ""
    pos_last = 0
    df_mutation$Chr_original = df_mutation$Chr
    if (!file.exists(bam_file_slim)) {
      for (mut_no in seq_len(dim(df_mutation)[1])) {
        print(paste(mut_no, "/", dim(df_mutation)[1]))
        if (mut_no == 1 & mut_no != dim(df_mutation)[1]) {
          if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
              df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
            syscom = paste("samtools view -h ",
                           bam_file,
                           " ",
                           df_mutation$Chr_original[mut_no],
                           ":",
                           max(1, df_mutation$Pos[mut_no] - 200),
                           "-",
                           df_mutation$Pos[mut_no] + 200,
                           " > ",
                           bam_file_slim,
                           sep="")
            system(syscom)
            continuous = FALSE
            sep_new = FALSE
          } else {
            continuous = TRUE
            pos_last = max(1, df_mutation$Pos[mut_no] - 200)
          }
        } else if (mut_no == 1 & mut_no == dim(df_mutation)[1]) {
          syscom = paste("samtools view -h ",
                         bam_file,
                         " ",
                         df_mutation$Chr[mut_no],
                         ":",
                         max(1, df_mutation$Pos[mut_no] - 200),
                         "-",
                         df_mutation$Pos[mut_no] + 200,
                         " > ",
                         bam_file_slim,
                         sep="")
          system(syscom)
        } else if (mut_no == dim(df_mutation)[1]) {
          if (sep_new) {
            syscom = paste("samtools view -h ",
                           bam_file,
                           " ",
                           df_mutation$Chr[mut_no],
                           ":",
                           pos_last,
                           "-",
                           df_mutation$Pos[mut_no] + 200,
                           " > ",
                           bam_file_slim,
                           sep="")
            system(syscom)
          } else if (continuous) {
            syscom = paste("samtools view -h ",
                           bam_file,
                           " ",
                           df_mutation$Chr_original[mut_no],
                           ":",
                           pos_last,
                           "-",
                           df_mutation$Pos[mut_no] + 200,
                           " > ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("samtools merge -fucp ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           " ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("mv ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           sep="")
            system(syscom)
          } else {
            syscom = paste("samtools view -h ",
                           bam_file,
                           " ",
                           df_mutation$Chr_original[mut_no],
                           ":",
                           max(1, df_mutation$Pos[mut_no] - 200),
                           "-",
                           df_mutation$Pos[mut_no] + 200,
                           " > ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("samtools merge -fucp ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           " ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("mv ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           sep="")
            system(syscom)
          }
        } else {
          if (sep_new) {
            if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
                df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
              syscom = paste("samtools view -h ",
                             bam_file,
                             " ",
                             df_mutation$Chr_original[mut_no],
                             ":",
                             pos_last,
                             "-",
                             df_mutation$Pos[mut_no] + 200,
                             " > ",
                             bam_file_slim,
                             sep="")
              system(syscom)
              continuous = FALSE
              sep_new = FALSE
            }
          } else if (continuous) {
            if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
                df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
              syscom = paste("samtools view -h ",
                             bam_file,
                             " ",
                             df_mutation$Chr_original[mut_no],
                             ":",
                             pos_last,
                             "-",
                             df_mutation$Pos[mut_no] + 200,
                             " > ",
                             bam_file_tmp1,
                             sep="")
              system(syscom)
              syscom = paste("samtools merge -fucp ",
                             bam_file_tmp2,
                             " ",
                             bam_file_slim,
                             " ",
                             bam_file_tmp1,
                             sep="")
              system(syscom)
              syscom = paste("mv ",
                             bam_file_tmp2,
                             " ",
                             bam_file_slim,
                             sep="")
              system(syscom)
              continuous = FALSE
            }
          } else {
            if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
                df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
              syscom = paste("samtools view -h ",
                             bam_file,
                             " ",
                             df_mutation$Chr_original[mut_no],
                             ":",
                             max(1, df_mutation$Pos[mut_no] - 200),
                             "-",
                             df_mutation$Pos[mut_no] + 200,
                             " > ",
                             bam_file_tmp1,
                             sep="")
              system(syscom)
              syscom = paste("samtools merge -fucp ",
                             bam_file_tmp2,
                             " ",
                             bam_file_slim,
                             " ",
                             bam_file_tmp1,
                             sep="")
              system(syscom)
              syscom = paste("mv ",
                             bam_file_tmp2,
                             " ",
                             bam_file_slim,
                             sep="")
              system(syscom)
            } else {
              continuous = TRUE
              pos_last = max(1, df_mutation$Pos[mut_no] - 200)
            }
          }
        }
      }
      syscom = paste("samtools view -bS ",
                     bam_file_slim,
                     " > ",
                     bam_file_tmp2,
                     sep="")
      system(syscom)
      syscom = paste("samtools sort -@ 4 -o ",
                     bam_file_slim,
                     " ",
                     bam_file_tmp2,
                     sep="")
      system(syscom)
      syscom = paste("rm ",
                     bam_file_tmp1,
                     sep="")
      system(syscom)
      syscom = paste("rm ",
                     bam_file_tmp2,
                     sep="")
      system(syscom)
      syscom = paste("samtools index ",
                     bam_file_slim,
                     sep="")
      system(syscom)
      print(paste("Slimmed BAM files were saved as ", bam_file_slim, sep=""))
    }
    
    bam_file = bam_file_slim
    fun_load_bam(bam_file) # df_bam
    if (list_exist) {
      fun_load_id(read_list) # df_mut_call
    }
    
    # analysis
    result = fun_read_check(short_homology_search_length = 4)
    msec = rbind(msec, result[[1]])
    homology_search = rbind(homology_search, result[[2]])
    mut_depth = list(rbind(mut_depth[[1]], result[[3]][[1]]),
                     rbind(mut_depth[[2]], result[[3]][[2]]),
                     rbind(mut_depth[[3]], result[[3]][[3]])) 
  }
  # search homologous sequences
  msec = fun_homology(msec,
                      homology_search,
                      min_homology_search = 40)
  
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
  list_exist <- TRUE
  
  # load genomic sequence
  fun_load_genome(organism) # ref_genome
  fun_load_chr_no(organism) # chr_no

  # load mutation information
  fun_load_mutation_gz(mutation_file) # df_mutation
  df_mutation = df_mutation[order(df_mutation$Chr, df_mutation$Pos),]
  sep_new = TRUE
  continuous = FALSE
  chr_last = ""
  pos_last = 0
  bam_file_tmp1 = paste(bam_file, ".tmp1", sep="")
  bam_file_tmp2 = paste(bam_file, ".tmp2", sep="")
  bam_file_slim = paste(bam_file, ".SLIM", sep="")
  if (!file.exists(bam_file_slim)) {
    for (mut_no in seq_len(dim(df_mutation)[1])) {
      print(paste(mut_no, "/", dim(df_mutation)[1]))
      if (mut_no == 1 & mut_no != dim(df_mutation)[1]) {
        if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
            df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
          syscom = paste("samtools view -h ",
                         bam_file,
                         " ",
                         df_mutation$Chr_original[mut_no],
                         ":",
                         max(1, df_mutation$Pos[mut_no] - 200),
                         "-",
                         df_mutation$Pos[mut_no] + 200,
                         " > ",
                         bam_file_slim,
                         sep="")
          system(syscom)
          continuous = FALSE
          sep_new = FALSE
        } else {
          continuous = TRUE
          pos_last = max(1, df_mutation$Pos[mut_no] - 200)
        }
      } else if (mut_no == 1 & mut_no == dim(df_mutation)[1]) {
        syscom = paste("samtools view -h ",
                       bam_file,
                       " ",
                       df_mutation$Chr_original[mut_no],
                       ":",
                       max(1, df_mutation$Pos[mut_no] - 200),
                       "-",
                       df_mutation$Pos[mut_no] + 200,
                       " > ",
                       bam_file_slim,
                       sep="")
        system(syscom)
      } else if (mut_no == dim(df_mutation)[1]) {
        if (sep_new) {
          syscom = paste("samtools view -h ",
                         bam_file,
                         " ",
                         df_mutation$Chr_original[mut_no],
                         ":",
                         pos_last,
                         "-",
                         df_mutation$Pos[mut_no] + 200,
                         " > ",
                         bam_file_slim,
                         sep="")
          system(syscom)
        } else if (continuous) {
          syscom = paste("samtools view -h ",
                         bam_file,
                         " ",
                         df_mutation$Chr_original[mut_no],
                         ":",
                         pos_last,
                         "-",
                         df_mutation$Pos[mut_no] + 200,
                         " > ",
                         bam_file_tmp1,
                         sep="")
          system(syscom)
          syscom = paste("samtools merge -fucp ",
                         bam_file_tmp2,
                         " ",
                         bam_file_slim,
                         " ",
                         bam_file_tmp1,
                         sep="")
          system(syscom)
          syscom = paste("mv ",
                         bam_file_tmp2,
                         " ",
                         bam_file_slim,
                         sep="")
          system(syscom)
        } else {
          syscom = paste("samtools view -h ",
                         bam_file,
                         " ",
                         df_mutation$Chr_original[mut_no],
                         ":",
                         max(1, df_mutation$Pos[mut_no] - 200),
                         "-",
                         df_mutation$Pos[mut_no] + 200,
                         " > ",
                         bam_file_tmp1,
                         sep="")
          system(syscom)
          syscom = paste("samtools merge -fucp ",
                         bam_file_tmp2,
                         " ",
                         bam_file_slim,
                         " ",
                         bam_file_tmp1,
                         sep="")
          system(syscom)
          syscom = paste("mv ",
                         bam_file_tmp2,
                         " ",
                         bam_file_slim,
                         sep="")
          system(syscom)
        }
      } else {
        if (sep_new) {
          if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
              df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
            syscom = paste("samtools view -h ",
                           bam_file,
                           " ",
                           df_mutation$Chr_original[mut_no],
                           ":",
                           pos_last,
                           "-",
                           df_mutation$Pos[mut_no] + 200,
                           " > ",
                           bam_file_slim,
                           sep="")
            system(syscom)
            continuous = FALSE
            sep_new = FALSE
          }
        } else if (continuous) {
          if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
              df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
            syscom = paste("samtools view -h ",
                           bam_file,
                           " ",
                           df_mutation$Chr_original[mut_no],
                           ":",
                           pos_last,
                           "-",
                           df_mutation$Pos[mut_no] + 200,
                           " > ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("samtools merge -fucp ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           " ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("mv ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           sep="")
            system(syscom)
            continuous = FALSE
          }
        } else {
          if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
              df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
            syscom = paste("samtools view -h ",
                           bam_file,
                           " ",
                           df_mutation$Chr_original[mut_no],
                           ":",
                           max(1, df_mutation$Pos[mut_no] - 200),
                           "-",
                           df_mutation$Pos[mut_no] + 200,
                           " > ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("samtools merge -fucp ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           " ",
                           bam_file_tmp1,
                           sep="")
            system(syscom)
            syscom = paste("mv ",
                           bam_file_tmp2,
                           " ",
                           bam_file_slim,
                           sep="")
            system(syscom)
          } else {
            continuous = TRUE
            pos_last = max(1, df_mutation$Pos[mut_no] - 200)
          }
        }
      }
    }
    syscom = paste("samtools view -bS ",
                   bam_file_slim,
                   " > ",
                   bam_file_tmp2,
                   sep="")
    system(syscom)
    syscom = paste("samtools sort -@ 4 -o ",
                   bam_file_slim,
                   " ",
                   bam_file_tmp2,
                   sep="")
    system(syscom)
    syscom = paste("rm ",
                   bam_file_tmp1,
                   sep="")
    system(syscom)
    syscom = paste("rm ",
                   bam_file_tmp2,
                   sep="")
    system(syscom)
    syscom = paste("samtools index ",
                   bam_file_slim,
                   sep="")
    system(syscom)
    print(paste("Slimmed BAM files were saved as ", bam_file_slim, sep=""))
  }
  
  bam_file = bam_file_slim
  fun_load_bam(bam_file) # df_bam
  fun_load_id(read_list) # df_mut_call

  # analysis
  result = fun_read_check(short_homology_search_length = 4)
  msec = result[[1]]
  homology_search = result[[2]]
  mut_depth = result[[3]]
  
  # search homologous sequences
  msec = fun_homology(msec,
                      homology_search,
                      min_homology_search = 40)
  
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

