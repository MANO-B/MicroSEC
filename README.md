![MicroSEC logo](MicroSEC_logo.png)

# Version history
2.1.3: Example files renewed.  
2.1.1: Bug fixed.  
2.1.0: Compatible with large indels in very short leads.  
2.0.0: Processing time reduced to 30% due to changes in search algorithm.  
  
# Docker and Apptainer
MicroSEC docker file can be downloadable via Docker-hub.  
https://hub.docker.com/r/ikegamitky/microsec/tags
```
docker pull ikegamitky/microsec:v2.1.3
```

Apptainer cantainer can be built via Docker-hub.
```
apptainer pull docker://ikegamitky/microsec:v2.1.3
```
  
Apptainer container can be built with a definition file (takes 30 min).  
```
wget https://raw.githubusercontent.com/MANO-B/MicroSEC/main/MicroSEC.def
sudo apptainer build MicroSEC.sif MicroSEC.def
```

  
# Analysis script file
MicroSEC has been improved to dramatically reduce memory usage.  
The speed of analysis is also improved by deleting parts of the BAM file that are not relevant to the mutations prior to analysis by MicroSEC.  
Samtools is now mandatory.  
Please download and use the new version of MicroSEC.R.  
```
wget https://raw.githubusercontent.com/MANO-B/MicroSEC/main/MicroSEC.R
```

# Known limitation
PIK3CA E545A (chr3:179218304A>C, NM_006218.4:c.1634A>C) pathogenic mutation might
be called as an artifact by MicroSEC, which may be a false positive error. 
A PIK3CA pseudogene (LOC100422375) in chromosome 22 harbors a base substitution 
which was identical with the mutation in PIK3CA c.1634A>C.  
Please check the reads manually with IGV.  

# MicroSEC pipeline for FFPE artifacts
This pipeline is designed for filtering sequence errors found in formalin-fixed and 
paraffin-embedded (FFPE) samples.  
This repository contains all the codes to regenerate results from our 
paper:  
"MicroSEC: Sequence error filtering pipeline for formalin-fixed and 
paraffin-embedded samples"  

## Supplementary Code
> M. Ikegami et al., "MicroSEC filters sequence errors for formalin-fixed and paraffin-embedded samples", Commun Biol 4, 1396 (2021). https://doi.org/10.1038/s42003-021-02930-4

## Contents

- [Overview](#overview)
- [Filtering detail](#filtering-detail)
- [System Requirements](#system-requirements)
- [Instructions for Use](#instructions-for-use)

## Overview

This pipeline is designed for filtering mutations found in formalin-fixed and paraffin-embedded (FFPE) samples.  
The MicroSEC filter utilizes a statistical analysis, and the results for mutations with less than 10 supporting reads are not reliable.  
Two files are necessary for the analysis: mutation information file, BAM file.  
A mutation supporting read ID information file is desirable but not necessary.  
Prepare a sample information tsv file.  
### File 1: mutation information tsv file (required)  
This tsv file should contain at least these contents (any number of other columns are allowed):  
```
Sample              Mut_type   Chr  Pos               Ref Alt SimpleRepeat_TRF Neighborhood_sequence  
SL_1010-N6-B  1-snv         chr1 108130741   C    T   N                               CTACCTGGAGAATGGGCCCATGTGTCCAGGTAGCAGTAAGC
```  
- Notation:  
    - Mut_type: [altered bases]-[mutation type (snv, ins, or del)]. 1-snv, 3-snv, 3-del, 1-ins, etc.  
    - Chr: chr1-22, chrX, and chrY for hg19 and hg38. 1-20, X, and Y for mm10.  
           MicroSEC handle both "chr1" and "1" format chromosomes.  
    - Pos/Ref/Alt: There are some difference from HGVS nomenclature as follows:  
          - 131C>T -> Pos/Ref/Alt 131/C/T.  
          - 514_515insC and 514A -> Pos/Ref/Alt 514/A/AC  
          - 244delC and 243A -> Pos/Ref/Alt 243/AC/A  
    - SimpleRepeat_TRF: Whether the mutation locates at a simple repeat sequence or not (Y or N).  
    - Neighborhood_sequence: [5'-20 bases] + [Alt sequence] + [3'-20 bases].  
    - Sample, Chr, Pos, Ref, and Alt should be set exactly.  
    - If you do not know Mut_type, SimpleRepeat_TRF, or Neighborhood_sequence, enter "-".
### File 2: BAM file (required)  
This file should contain at least these contents (always included in standard BAM files):  
- QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, ISIZE, and QUAL.  

I only have experience processing BAM files of up to 20 gigabytes.  
It is recommended that huge BAM files be split up for processing.  
  
For example,  
```
samtools view -b sort.bam chr1:100-1000 > output_chr1.bam
```

It is better to split the file by chromosome.  
If the file size is still too large, split it further at a distance from the mutation.  
Deleting regions where there are no mutations will lighten the process.  
Please download and use MicroSEC.R script file.  
  
### File 3: sample information tsv file (required)  
From seven to ten columns are necessary. Three optional columns can be omitted.  
The file contains no header.
```
[sample name] [mutation information tsv file] [BAM file] [read length] [adapter sequence read 1] [optional: adapter sequence read 2] [sample type: Human or Mouse] [panel name] [optional: reference genome fastq file] [optional: simple repeat region bed file]  
PC9	./source/CCLE.tsv	./source/Cell_line/PC9.bam 127	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT Human TOP
A375 ./source/CCLE.tsv.gz	./source/Cell_line/A375.bam	127	 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT Hg38 TOP ./reference/hg38.fa ./reference/simpleRepeat.bed.gz
```
  
- Reference genome: Human (hg38), Mouse (mm10), hg19, hg38, or mm10.
  Support for other organisms is easy, and will be done as needed upon request.  

### File 4: Reference genome  
Reference fastq file is necessary when you using CRAM files.  
  
### File 5: simple repeat region bed file (optional, but mandatory to detect simple repeat derived artifacts)  
simpleRepeat.bed.gz file is necessary when you check simple repeat derived errors.  
    
## Filtering detail
This pipeline contains 8 filtering processes.  

- Filter 1  : Shorter-supporting lengths distribute too short to occur (1-1 and 1-2).  
        - Filter 1-1: P-values are less than the threshold_p(default: 10^(-6)).  
        - Filter 1-2: The shorter-supporting lengths distributed over less than 75% of the read length.  
- Filter 2  : Hairpin-structure induced error detection (2-1 and 2-2).  
        - Filter 2-1: Palindromic sequences exist within 200 bases.  
        - Filter 2-2: >=50% mutation-supporting reads contains a reverse complementary sequence of the opposite strand consisting >= 15 bases.  
- Filter 3  : 3'-/5'-supporting lengths are too densely distributed to occur (3-1 and 3-2).  
        - Filter 3-1: P-values are less than the threshold_p(default: 10^(-6)).  
        - Filter 3-2: The distributions of 3'-/5'-supporting lengths are within 75% of the read length.  
- Filter 4  : >=15% mutations were called by chimeric reads comprising two distant regions.  
- Filter 5  : >=50% mutations were called by soft-clipped reads.  
- Filter 6  : Mutations locating at simple repeat sequences.  
- Filter 7  : Indel mutations locating at a >=15 homopolymer.  
- Filter 8  : >=10% of bases are low quality (Quality score <18) in the mutation supporting reads.  

Filter 1, 2, 3, and 4 detect possible FFPE artifacts.  
Filter 5 may also be FFPE artifacts or mapping errors.  
Filter 6, 7, and 8 detect frequent errors caused by the next generation sequencing platform.  
Supporting lengths are adjusted considering small repeat sequences around the mutations.  
  
Results are saved in a tsv file.  

github url: https://github.com/MANO-B/MicroSEC

## System Requirements

### Hardware Requirements

The scripts requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 4 GB of RAM (depending on the size of BAM file and the number of mutations). For optimal performance, we recommend a computer with the following specs:

RAM: 4+ GB  
CPU: 2+ cores, 2.6+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, M1 Macbook air) and internet of speed 40 Mbps.

### Software Requirements

### Samtools

Samtools is used for pre-processing to remove reads that are not related to mutations. Version 1.12 is what I am using, but I think older versions will work if they support multi-core processing.  

### R language

This script files runs on `R` for Windows, Mac, or Linux, which requires the R version 3.5.0 or later.
If you use version 3.4 or lower of R, you will have some difficulty installing the packages, but it is not impossible.  

### Package dependencies

Users should install the following packages prior to use the scripts, from an `R` terminal:

```
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
install.packages(c('stringr', 'dplyr', 'remotes'), dependencies = TRUE)
BiocManager::install(c("Rsamtools", "Biostrings", "GenomicAlignments", "GenomeInfoDb"), update=FALSE)

# install necessary genomes
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=FALSE)
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", update=FALSE)
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update=FALSE)
```

which will install in about 30 minutes on a recommended machine.

### Package Versions
The program does not use any of the functions specific to the following version of the packages, so there is no problem if you use the latest version of the package.  

```
> packageVersion("stringr")
[1] '1.4.0'
> packageVersion("dplyr")
[1] '1.0.2'
> packageVersion("Biostrings")
[1] '2.54.0'
> packageVersion("BSgenome.Hsapiens.UCSC.hg38")
[1] '1.4.1'
> packageVersion("GenomicAlignments")
[1] '1.22.1'
> packageVersion("Rsamtools")
[1] '2.0.3'
> packageVersion("remotes")
[1] '2.5.0'
> packageVersion("GenomeInfoDb")
[1] '1.22.1'
```

## Instructions for Use
See also https://rdrr.io/cran/MicroSEC/
- How to install
```
# Stable version (v2.1.3) from github (recommended)
remotes::install_github("MANO-B/MicroSEC", ref = 'v2.1.3')

# Developmental version from github
remotes::install_github("MANO-B/MicroSEC")
```
- How to use in command line
```
# download only once
wget https://raw.githubusercontent.com/MANO-B/MicroSEC/main/MicroSEC.R
# run the script
Rscript MicroSEC.R [working/output directory] [sample information tsv file] [progress bar Y/N]
```  
- Example
```
Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/Sample_list.txt Y  
```  
- How to use in R Console
```
## See also the vignette file.

## Necessary packages
library(MicroSEC)

options(show.error.messages = FALSE, warn = -1)
wd <- "mnt/HDD8TB/MicroSEC" # set your working/output directory
sample_list <- "mnt/HDD8TB/MicroSEC/source/Sample_list.txt"
progress_bar <- "Y"
  
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
  bam_file_slim <- paste0(bam_file, ".SLIM")
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
  result <- fun_read_check(df_mutation = exampleMutation,
                           df_bam = exampleBAM,
                           ref_genome = ref_genome,
                           sample_name = sample_name,
                           read_length = read_length,
                           adapter_1 = adapter_1,
                           adapter_2 = adapter_2,
                           short_homology_search_length = 4,
                           min_homology_search = 40,
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
fun_save(msec, paste0(wd, "/", sample_info[1,1], ".tsv"))
```

- Output files  
  A tsv file is saved in the working/output directory.
```
[Sample name].tsv
```
  
- Sample name is set to the sample of interest in the Sample column of the mutation information file.
- Confirm the read length in the platform.
- Confirm the adapter sequence; Illumina sequencer ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" and "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT").
- If you input only one adapter sequence, the sequence will be used for read 1 and read 2. 
- If you want to know the progress visually, [progress bar Y/N] should be Y.
  
### Reproducibility

The source code is available in MicroSEC.R. 
Source data will be available at the Japanese 
Genotype-Phenotype Archive (http://trace.ddbj.nig.ac.jp/jga), 
which is hosted by the DNA Data Bank of Japan (the accession 
number is written in our paper). Each filtering process takes 
about 5–120 minutes per sample on a recommended machine, 
according to the depth and mutation amount of the sample.
