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
# This tsv or tsv.gz file should contain at least these columns:
#       Sample Mut_type   Chr       Pos Ref Alt SimpleRepeat_TRF                     Neighborhood_sequence
# SL_1010-N6-B    1-snv  chr1 108130741   C   T                N CTACCTGGAGAATGGGCCCATGTGTCCAGGTAGCAGTAAGC
# SimpleRepeat_TRF: Whether the mutation locates at a simple repeat sequence or not ("Y" or "N").
# Neighborhood_sequence: [5'-20 bases] + [Alt sequence] + [3'-20 bases]
# Sample, Mut_type, Chr, Pos, Ref, and Alt should be set exactly.  
# If you do not know the SimpleRepeat_TRF, Mut_type, or Neighborhood_sequence, enter "-". Automatically detected.
#
# File 2: BAM file (mandatory)
#
# File 3: sample information tsv file  (mandatory, if multiple samples are processed in a batch)
# Seven to ten columns are necessary (without column names).
# Optional columns can be deleted if they are not applicable.
# [sample name] [mutation information tsv file] [BAM file] [read length] [adapter sequence read 1] [optional: adapter sequence read 2] [sample type: Human or Mouse] [panel name] [optional: reference genome fastq file] [optional: simple repeat region bed file]  
# PC9	./source/CCLE.tsv	./source/Cell_line/PC9.bam 127	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT Human TOP
# A375 ./source/CCLE.tsv.gz	./source/Cell_line/A375.bam	127	 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT Hg38 TOP ./reference/hg38.fa ./reference/simpleRepeat.bed.gz
#
# File 4: Reference genome: Human (hg38 or hg19) or Mouse (mm10) (optional, but mandatory with cram files)
#
# File 5: simple repeat region bed file (optional, but mandatory to detect simple repeat derived artifacts)  
#
# This pipeline contains 8 filtering processes.
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
# Filter 5  : >=50% mutations were called by soft-clipped reads.
# Filter 6  : Mutations locating at simple repeat sequences.
# Filter 7  : Mutations locating at a >=15 homopolymer.
# Filter 8  : >=10% low quality bases (Quality score <18) in the mutation supporting reads.
#
# Supporting lengths are adjusted considering small repeat sequences around the mutations.
#
# How to use
# Rscript MicroSEC.R [working/output directory] [sample information tsv file] [progress bar Y/N]
# 
# Example
# Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/Sample_list_test.txt N
# Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source//sample_info_test.tsv.gz Y
#
# If you want to know the progress visually, [progress bar Y/N] should be Y.
#
# Results are saved in a tsv file.

# load necessary packages
library(MicroSEC)
library(dplyr)
library(openxlsx)
library(readr)
library(stringr)
library(Rsamtools)
library(BiocGenerics)
library(lubridate)
library(Biostrings)
library(gtools)

# set arguments
args <- commandArgs(trailingOnly = T)

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
  sample_name <- sample_info[sample, 1]
  mutation_file <- sample_info[sample, 2]
  bam_file <- sample_info[sample, 3]
  read_length <- as.integer(sample_info[sample, 4])
  adapter_1 <- sample_info[sample, 5]
  if (sample_info[sample, 6] %in%
      c("Human", "Mouse", "hg19", "hg38", "mm10")) {
    adapter_2 <- adapter_1
    organism <- sample_info[sample, 6]
    panel <- sample_info[sample, 7]
    if (dim(sample_info)[2] == 8) {
      reference_genome <- sample_info[sample, 8]
    }
    if (dim(sample_info)[2] == 9) {
      reference_genome <- sample_info[sample, 8]
      simple_repeat_list <- sample_info[sample, 9]
    }
  } else{
    adapter_2 <- sample_info[sample, 6]
    organism <- sample_info[sample, 7]
    panel <- sample_info[sample, 8]
    if (dim(sample_info)[2] == 9) {
      reference_genome <- sample_info[sample, 9]
    }
    if (dim(sample_info)[2] == 10) {
      reference_genome <- sample_info[sample, 9]
      simple_repeat_list <- sample_info[sample, 10]
    }
  }
  bam_file_slim <- paste0(bam_file, ".SLIM.bam")
  bam_file_tmp = paste0(bam_file, ".tmp.bam")
  
  # load genomic sequence
  ref_genome <- fun_load_genome(organism)
  chr_no <- fun_load_chr_no(organism)
  if (ref_genome@user_seqnames[[1]] == "chr1") {
    chromosomes <- paste0("chr", c(seq_len(chr_no - 2),"X", "Y"))
  }
  if (ref_genome@user_seqnames[[1]] == "1") {
    chromosomes <- paste0("", c(seq_len(chr_no - 2),"X", "Y"))
  }
  
  # load mutation information
  df_mutation <- fun_load_mutation(mutation_file, sample_name, ref_genome)
  
  if (tools::file_ext(bam_file) == "bam") {
    bam_file_bai <- paste0(bam_file, ".bai")
  } else if (tools::file_ext(bam_file) == "cram") {
    bam_file_bai <- paste0(bam_file, ".crai")
  }
  if (!file.exists(bam_file_bai) & !file.exists(bam_file_slim)) {
    print("Sorting a BAM/CRAM file...")
    if (tools::file_ext(bam_file) == "bam") {
      bam_file_sort <- paste0(bam_file, "_sort.bam")
      syscom <- paste0("samtools sort -@ 4 -o ",
                       bam_file_sort, " ",
                       bam_file)
    } else if (tools::file_ext(bam_file) == "cram") {
      bam_file_sort <- paste0(bam_file, "_sort.cram")
      syscom <- paste0("samtools sort -@ 4 -O cram -o ",
                       bam_file_sort," ",
                       bam_file)
    }        
    system(syscom)
    syscom <- paste0("samtools index ",
                    bam_file_sort)
    system(syscom)
    bam_file <- bam_file_sort
  }

  if (!file.exists(paste0(bam_file_slim, ".bai"))) {
    print("Trimming a BAM/CRAM file...")
    df_mutation_save <- df_mutation
    download_region <- data.frame(matrix(rep(NA,3),nrow=1))[numeric(0),]
    colnames(download_region) <- c("chrom", "chromStart", "chromEnd")
    print(paste0(sample_name, ": ", dim(df_mutation_save)[1], " mutations"))
    for (i in chromosomes) {
      df_mutation <- df_mutation_save[df_mutation_save$Chr == i,]
      continuous = FALSE
      for (mut_no in seq_len(dim(df_mutation)[1])) {
        if (mut_no == 1 & mut_no != dim(df_mutation)[1]) {
          if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
             df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
            download_region <- rbind(
              download_region,
              c(df_mutation$Chr[mut_no],
                max(1, df_mutation$Pos[mut_no] - 200) - 1,
                df_mutation$Pos[mut_no] + 200))
            colnames(download_region) <-
              c("chrom", "chromStart", "chromEnd")
            continuous <- FALSE
          } else {
            continuous <- TRUE
            pos_last <- max(1, df_mutation$Pos[mut_no] - 200)
          }
        } else if (mut_no == 1 & mut_no == dim(df_mutation)[1]) {
          download_region <- rbind(
            download_region,
            c(df_mutation$Chr[mut_no],
              max(1, df_mutation$Pos[mut_no] - 200) - 1,
              df_mutation$Pos[mut_no] + 200))
          colnames(download_region) <-
            c("chrom", "chromStart", "chromEnd")
        } else if (mut_no == dim(df_mutation)[1]) {
          if (continuous) {
            download_region <- rbind(
              download_region,
              c(df_mutation$Chr[mut_no],
                pos_last - 1,
                df_mutation$Pos[mut_no] + 200))
            colnames(download_region) <-
              c("chrom", "chromStart", "chromEnd")
          } else {
            download_region = rbind(
              download_region,
              c(df_mutation$Chr[mut_no],
                max(1, df_mutation$Pos[mut_no] - 200) - 1,
                df_mutation$Pos[mut_no] + 200))
            colnames(download_region) <-
              c("chrom", "chromStart", "chromEnd")
          }
        } else {
          if (continuous) {
            if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
                df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
              download_region = rbind(
                download_region,
                c(df_mutation$Chr[mut_no],
                  pos_last - 1,
                  df_mutation$Pos[mut_no] + 200))
              colnames(download_region) <-
                c("chrom", "chromStart", "chromEnd")
              continuous = FALSE
            }
          } else {
            if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
               df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
              download_region <- rbind(
                download_region,
                c(df_mutation$Chr[mut_no],
                  max(1, df_mutation$Pos[mut_no] - 200) - 1,
                  df_mutation$Pos[mut_no] + 200))
              colnames(download_region) <-
                c("chrom", "chromStart", "chromEnd")
            } else {
              continuous <- TRUE
              pos_last <- max(1, df_mutation$Pos[mut_no] - 200)
            }
          }
        }
      }
    }
    write_tsv(x = download_region,
              file = paste0(bam_file,".bed"),
              progress = F,
              col_names = F)
    rm(download_region)
    system_out <- 1
    if (tools::file_ext(bam_file) == "bam") {
      syscom <- paste0("samtools view -h --no-PG ",
                     bam_file, " -L ",
                     paste0(bam_file,".bed"), " > ",
                     bam_file_slim)
    } else if (tools::file_ext(bam_file) == "cram") {
      file_crai <- paste0(bam_file, ".crai")
      syscom <- paste0("samtools view -h --no-PG ",
                       bam_file, " -T ",
                       reference_genome, " ",
                       "-X ", file_crai, " -L ",
                       paste0(bam_file,".bed"), " > ",
                       bam_file_slim)
    }
    system_out = (system(syscom))
    if(system_out == 0){
      system_out <- 1
      syscom <- paste0("samtools sort -o ",
                       bam_file_tmp, " ",
                       bam_file_slim)
      print("Sorting BAM file...")
      system_out <- (system(syscom))
      if(system_out == 0){
        system_out <- 1
        syscom = paste0("samtools view -bS ",
                        bam_file_tmp, " > ",
                        bam_file_slim)
        print("Compressing BAM file...")
        system_out <- (system(syscom))
        if(system_out == 0){
          system_out <- 1
          syscom = paste0("samtools index ", bam_file_slim)
          print("Indexing BAM file...")
          system_out <- (system(syscom))
          if(system_out == 0){
            print(paste0("Slimmed BAM files were saved as ", bam_file_slim))
            #file.remove(bam_file)
            file.remove(bam_file_tmp)
          }
        }
      }
    }
    df_mutation <- df_mutation_save
    rm(df_mutation_save)
  }

  bam_file <- bam_file_slim
  df_bam <- fun_load_bam(bam_file)

  # analysis
  result <- fun_read_check(df_mutation = df_mutation,
                           df_bam = df_bam,
                           ref_genome = ref_genome,
                           sample_name = sample_name,
                           read_length = read_length,
                           adapter_1 = adapter_1,
                           adapter_2 = adapter_2,
                           short_homology_search_length = 4,
                           progress_bar = progress_bar)
  msec <- rbind(msec, result[[1]])
  homology_search <- rbind(homology_search, result[[2]])
  mut_depth <- list(rbind(mut_depth[[1]], result[[3]][[1]]),
                   rbind(mut_depth[[2]], result[[3]][[2]]),
                   rbind(mut_depth[[3]], result[[3]][[3]])) 
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
                    threshold_soft_clip_ratio = 0.50,
                    threshold_low_quality_rate = 0.1,
                    homopolymer_length = 15)

# save the results
fun_save(msec, paste0(wd, "/", sample_info[1,1], ".tsv.gz"))
