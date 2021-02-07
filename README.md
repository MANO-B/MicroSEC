![MicroSEC logo](MicroSEC_logo.png)

# MicroSEC pipeline for FFPE artifacts
This pipeline is designed for filtering sequence errors found in formalin-fixed and 
paraffin-embedded (FFPE) samples.  
This repository contains all the codes to regenerate results from our 
paper:  
"MicroSEC: Sequence error filtering pipeline for formalin-fixed and 
paraffin-embedded samples"  

## Supplementary Code
> M. Ikegami et al., "MicroSEC: Sequence error filtering pipeline for 
formalin-fixed and paraffin-embedded samples", in preparation.
DOI: 10.21203/rs.3.rs-153650/v1
https://www.researchsquare.com/article/rs-153650/v1

## Contents

- [Overview](#overview)
- [Filtering detail](#filtering-detail)
- [System Requirements](#system-requirements)
- [Instructions for Use](#instructions-for-use)

## Overview

This pipeline is designed for filtering mutations found in formalin-fixed and paraffin-embedded (FFPE) samples.  
The MicroSEC filter utilizes a statistical analysis, and the results for mutations with less than 10 supporting reads are not reliable.  
Two files are necessary for the analysis: mutation information file, BAM file  
A mutation supporting read ID information file is desirable but not necessary.  
Prepare a sample information tsv file.  
### File 1: mutation information file  
This excel file should contain at least these contents:  
```
Sample       Gene       HGVS.c  HGVS.p Mut_type Total_QV>=20   %Alt   Chr  Pos         Ref Alt SimpleRepeat_TRF Neighborhood_sequence                      Transition  
SL_1010-N6-B SLC25A24    _      _      1-snv    366            1.0929 chr1 108130741   C   T   N                CTACCTGGAGAATGGGCCCATGTGTCCAGGTAGCAGTAAGC  C>T_t
```  
- Notation:  
    - Mut_type: [altered bases]-[mutation type (snv, ins, or del)]. 1-snv, 3-snv, 3-del, 1-ins, etc.  
    - Total_QV>=20: The read number with total Q-value >=20.  
    - Chr: chr1-22, chrX, and chrY for hg38. 1-22, X, and Y for hg19.  
    - Pos/Ref/Alt: There are some difference from HGVS nomenclature as follows:  
          - 131C>T -> Pos/Ref/Alt 131/C/T.  
          - 514_515insC and 514A -> Pos/Ref/Alt 514/A/AC  
          - 244delC and 243A -> Pos/Ref/Alt 243/AC/A  
    - SimpleRepeat_TRF: Whether the mutation locates at a simple repeat sequence or not.  
    - Neighborhood_sequence: [5'-20 bases] + [Alt sequence] + [3'-20 bases].  
    - Transition: 1-snv mutation pattern with a 3'-base. C>T_t represents CT to TT mutation. C>T_g_FFPE represents the possible FFPE artifact.  
    - Sample, Mut_type, Chr, Pos, Ref, and Alt should be set exactly.  
    - Gene, HGVS.c, HGVS.p, Total_QV>=20, %Alt, SimpleRepeat_TRF, and Transition can be set to any values.  
    - If you do not know the Neighborhood_sequence, enter "-".
### File 2: BAM file  
This file should contain at least these contents:  
- QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, and QUAL  
  
### File 3: mutation supporting read ID information tsv file  
This file should contain at least these contents:  
```
Chr  Pos       Ref Alt        Mut_ID                                                                                                Mut
chr1 2561609   T   A          _;ID001-1:579185f,ID004-1:1873933f;ID006-1:1131647f,ID001-1:570086f,ID008-1:1953407r,ID002-2:749570r  .;A;N
chr6 346429    G   GACACACAC  _;ID005-2:545593f,ID006-1:1132212f,ID006-1:1132212r;ID004-1:991476r,ID004-2:516896f;ID002-1:670026f   .;.+ACACACAC;.-ACAC;N
```
- Mut_ID  
    _;AAA,BBB,CCC;DDD,EEE,FFF represents the ID list corresponding to the Mut column separated by semicolons.  
    Read ID is "QNAME in the BAM file" + "mapped strand information".  
    Strand information can be extracted from the FLAG column in the BAM file.  
    ID001-1:579185f represents that the read ID is ID001-1:579185 and the read was mapped to the forward strand.  
    ID008-1:1953407r represents that the read ID is ID008-1:1953407 and the read was mapped to the reverse strand.  
- Mut  
    ".;A;N" represents the mutation pattern: "reference base";"altered to A";"other alterations"  
    ".;.+ACACACAC;.-ACAC;N" reapresents the mutation pattern: "reference base";"insertion ACACACAC";"deletion ACAC";"other alterations"  
  
### File 4: sample information tsv file  
From six to eight columns are necessary. Two optional columns can be omitted.  
The file contains no header.
```
[sample name] [mutation information excel file] [BAM file]          [optional: read ID information directory] [read length] [adapter sequence read 1]         [optional: adapter sequence read 2] [reference genome]  
PC9           /mnt/source/CCLE.xlsx             /mnt/source/PC9.bam /mnt/source/PC9_Cell_line                 127           AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT   hg38  
```
  
- Reference genome: Human (hg38), Mouse (mm10), hg19, hg38, or mm10
  
## Filtering detail
This pipeline contains 8 filtering processes.  

- Filter 1  : Shorter-supporting lengths distribute too short to occur (1-1 and 1-2).  
        - Filter 1-1: P-values are less than the threshold_p(default: 10^(-6)).  
        - Filter 1-2: The longest shorter-supporting lengths is shorter than 40% of the read length.  
- Filter 2  : Hairpin-structure induced error detection (2-1 or 2-2).  
        - Filter 2-1: Palindromic sequences exist within 150 bases.  
        - Filter 2-2: >=50% mutation-supporting reads contains a reverse complementary sequence of the opposite strand consisting >= 15 bases.  
- Filter 3  : 3'-/5'-supporting lengths are too densely distributed to occur (3-1, 3-2, and 3-3).  
        - Filter 3-1: P-values are less than the threshold_p(default: 10^(-6)).  
        - Filter 3-2: The distributions of 3'-/5'-supporting lengths are shorter than 80% of the read length.  
        - Filter 3-3: <10% of bases are low quality (Quality score <18).  
- Filter 4  : >=15% mutations were called by chimeric reads comprising two distant regions.  
- Filter 5  : Mutations locating at simple repeat sequences.  
- Filter 6  : C>T_g false positive calls in FFPE samples.  
- Filter 7  : Indel mutations locating at a >=15 homopolymer.  
- Filter 8  : >=10% low quality bases in the mutation supporting reads.  

Filter 1, 2, 3, and 4 detect possible FFPE artifacts.  
Filter 5, 6, 7, and 8 detect frequent errors caused by next generation sequencers.  
Supporting lengths are adjusted considering small repeat sequences around the mutations.  
  
Results are saved in a excel file.  
The explatation of the results is written in detail in the second sheet of the excel file.  

github url: https://github.com/MANO-B/MicroSEC

## System Requirements

### Hardware Requirements

The scripts requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 32 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 128+ GB  
CPU: 4+ cores, 4.2+ GHz/core

The runtimes below are generated using a computer with the recommended specs (128 GB RAM, 4 cores@4.2 GHz) and internet of speed 100 Mbps.

### Software Requirements

### R language

This script files runs on `R` for Windows, Mac, or Linux, which requires the R version 3.4.0 or later.


### Package dependencies

Users should install the following packages prior to use the scripts, from an `R` terminal:

```
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
install.packages(c('tidyr', 'openxlsx', 'data.table', 'R.utils', 'stringr', 'magrittr', 'dplyr', 'gtools', 'devtools'), dependencies = TRUE)
BiocManager::install(c("Rsamtools", "Biostrings", "GenomicAlignments", "GenomeInfoDb"), update=FALSE)

# install necessary genomes
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=FALSE)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", update=FALSE)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update=FALSE)
```

which will install in about 30 minutes on a recommended machine.

### Package Versions

All packages are in their latest versions as they appear on `CRAN` on Oct. 31, 2020. Users can check [CRAN snapshot](https://mran.microsoft.com/timemachine/) for details. The versions of software are, specifically:

```
> packageVersion("tidyr")
[1] '1.1.2'
> packageVersion("openxlsx")
[1] '4.2.2'
> packageVersion("data.table")
[1] '1.13.2'
> packageVersion("R.utils")
[1] '2.10.1'
> packageVersion("stringr")
[1] '1.4.0'
> packageVersion("magrittr")
[1] '1.5'
> packageVersion("dplyr")
[1] '1.0.2'
> packageVersion("Biostrings")
[1] '2.54.0'
> packageVersion("BSgenome.Hsapiens.UCSC.hg38")
[1] '1.4.1'
> packageVersion("BSgenome.Hsapiens.UCSC.hg19")
[1] '1.4.0'
> packageVersion("BSgenome.Mmusculus.UCSC.mm10")
[1] '1.4.0'
> packageVersion("GenomicAlignments")
[1] '1.22.1'
> packageVersion("Rsamtools")
[1] '2.0.3'
> packageVersion("gtools")
[1] '3.8.2'
> packageVersion("devtools")
[1] '2.3.2'
> packageVersion("GenomeInfoDb")
[1] '1.22.1'
```

## Instructions for Use
See also https://rdrr.io/cran/MicroSEC/
- How to install
```
# Stable version from CRAN
install.packages("MicroSEC")

# Developmental stable version from github (recommended)
devtools::install_github("MANO-B/MicroSEC", upgrade="never", ref="v1.2.1")

# Developmental unstable version from github (not recommended)
devtools::install_github("MANO-B/MicroSEC", upgrade="never")
```
- How to use in command line
```
# download only once
wget https://raw.githubusercontent.com/MANO-B/MicroSEC/main/MicroSEC.R

# Style 1
Rscript MicroSEC.R [working/output directory] [sample information tsv file] [progress bar Y/N]
 
# Style 2 (for pipeline use)
Rscript MicroSEC.R [output filename] [sample name] [mutation information tsv.gz] [BAM file] [read ID information directory] [read length] [adapter sequence read 1] [adapter sequence read 2] [sample type: Human or Mouse]
```  
- Example
```
Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/Sample_list.txt Y  
Rscript MicroSEC.R /mnt/result/post_filter/SAMPLE.gz SAMPLE_NAME /mnt/result/mutation/SAMPLE.gz /mnt/result/BAM/SAMPLE.bam /mnt/result/ID 150 	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT hg38
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

  # load genomic sequence
  fun_load_genome(organism)
  fun_load_chr_no(organism)

  # load mutation information
  fun_load_mutation(mutation_file, sample_name)
  fun_load_bam(bam_file)
  fun_load_id(read_list)

  # analysis
  result <- fun_read_check(short_homology_search_length = 4)
  msec <- rbind(msec, result[[1]])
  homology_search <- rbind(homology_search, result[[2]])
  mut_depth <- rbind(mut_depth, result[[3]])
}
# search homologous sequences
msec <- fun_homology(msec,
                    homology_search,
                    min_homology_search = 40)

# statistical analysis
msec <- fun_summary(msec)
msec <- fun_analysis(msec,
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
```

- Output files  
  one excel file and one tsv file are saved in the working/output directory.
```
MicroSEC-result_[Sample name]_[Date].xlsx
MicroSEC_[Sample name].tsv
```
  like following:  
```
MicroSEC-result_PC9_2020-11-05.xlsx
MicroSEC_PC9.tsv
```
  
- Sample name is set to the sample of interest in the Sample column of the mutation information file.
- Confirm the read length in the platform.
- Confirm the adapter sequence; Todai Onco Panel ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" and "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT").
- If you input only one adapter sequence, the sequence will be used for read 1 and read 2. 
- If you want to know the progress visually, [progress bar Y/N] should be Y.
  
  
### Conversion of mutation information files
MicroSEC utilizes a mutation information file. The format of the file is written differently from common HGVS nomenclature.  
A conversion function is prepared.
```
# Conversion example
# Set column names correctly.
#
# Before (hg19)
# Sample Hugo_Symbol    Chr   Start_position   End_position Variant_Type Reference Tumor_Seq Protein_Change
# LS411N BRAF           chr7     140453136        140453136    SNP          A         T         p.V600E
# LS411N CELSR1         chr22    46931227         46931227     DEL          C         -         p.G615fs
# LS411N APC            chr5     112175951        112175952    INS          -         A         p.E1554fs
#
# After
# Sample Gene     Chr   Pos          Mut_type  Ref   Alt HGVS.p       Neighborhood_sequence
# LS411N BRAF     chr7     140453136    1-snv     A     T   p.V600E      ACCCACTCCATCGAGATTTCTCTGTAGCTAGACCAAAATCA
# LS411N CELSR1   chr22    46931226     1-del     GC    G   p.G615fs     TCTTAGGCCCAGCGCTGCCGCCCCCAGAAAGGTGGAGGCC
# LS411N APC      chr5     112175951    1-ins     G     GA  p.E1554fs    AAAACCAAGAGAAAGAGGCAGAAAAAAACTATTGATTCTGAA    

library(openxlsx)
library(MicroSEC)

df_mutation = fun_convert(
  system.file("extdata", "convert_test.xlsx", package = "MicroSEC"),
  "hg19"
)

write.xlsx(df_mutation, "/mnt/HDD8TB/MicroSEC/mutation_modified.xlsx")
```

### Try to apply MicroSEC to CCLE WES data
CCLE (Cancer Cell Line Encyclopedia) supplies whole exome sequencing 
data and somatic mutation call lists of lots of cell lines.  
This is an open source dataset.  
Although these data are obtained from not FFPE samples but fresh samples, 
some possible MICR-originating sequence artifacts are detected.

- Download 24 BAM file from the NCBI SRA website.
  Choose "Output this run in BAM format to File".  
  The files will be downloaded as SRR8618961_1.bam - SRR8618961_Y.bam.  
  https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8618961  
- Merge all BAM files with samtools
```
a="SRR8618961"
samtools merge -@ 4 ${a}_all.bam ${a}_1.bam ${a}_2.bam ${a}_3.bam ${a}_4.bam ${a}_5.bam ${a}_6.bam ${a}_7.bam ${a}_8.bam ${a}_9.bam ${a}_10.bam ${a}_11.bam ${a}_12.bam ${a}_13.bam ${a}_14.bam ${a}_15.bam ${a}_16.bam ${a}_17.bam ${a}_18.bam ${a}_19.bam ${a}_20.bam ${a}_21.bam ${a}_22.bam ${a}_X.bam ${a}_Y.bam
samtools sort -@ 4 ${a}_all.bam > ${a}_sort.bam
samtools rmdup ${a}_sort.bam ${a}.bam 
samtools index -@ 4 ${a}.bam
```

- Download a mutation call list from CCLE website.
  "CCLE_DepMap_18q3_maf_20180718.txt"  
  https://portals.broadinstitute.org/ccle  
- Convert the mutation list with fun_convert function.
```
library(openxlsx)
library(MicroSEC)
library(GenomeInfoDb)

df_CCLE <- utils::read.csv("CCLE_DepMap_18q3_maf_20180718.txt",
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          sep = "\t")
df_CCLE <- df_CCLE %>% dplyr::filter(Tumor_Sample_Barcode == "NCIH2286_LUNG")
df_CCLE$Sample <- df_CCLE$Tumor_Sample_Barcode
df_CCLE$Chr <- df_CCLE$Chromosome
df_CCLE$Reference <- df_CCLE$Reference_Allele
df_CCLE$Tumor_Seq <- df_CCLE$Tumor_Seq_Allele1
seqlevelsStyle(df_CCLE$Chr) <- "UCSC"
df_CCLE <- df_CCLE %>% dplyr::filter(Chr != "M")
write.xlsx(df_CCLE, "/mnt/HDD8TB/MicroSEC/source/NCIH2286_LUNG_Mutations.xlsx")

df_mutation <- fun_convert(
  "NCIH2286_LUNG_Mutations.xlsx",
  "hg19"
)

write.xlsx(df_mutation, "/mnt/HDD8TB/MicroSEC/source/NCIH2286_LUNG_modified.xlsx")
```
- Perform MicroSEC filtering
```
msec <- NULL
homology_search <- NULL
mut_depth <- NULL
progress_bar <- "Y"

sample_name <- "NCIH2286_LUNG"
mutation_file <- "/mnt/HDD8TB/MicroSEC/source/NCIH2286_LUNG_modified.xlsx"
bam_file <- "SRR8618961.bam"
list_exist <- FALSE
read_length <- 76
adapter_1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
organism <- "hg19"

# load genomic sequence
fun_load_genome(organism) # ref_genome
fun_load_chr_no(organism) # chr_no

fun_load_mutation(mutation_file, sample_name) # df_mutation
fun_load_bam(bam_file) # df_bam

result <- fun_read_check(short_homology_search_length = 4)
msec <- rbind(msec, result[[1]])
homology_search <- rbind(homology_search, result[[2]])
mut_depth <- rbind(mut_depth, result[[3]])
msec <- fun_homology(msec,
                    homology_search,
                    min_homology_search = 40)

# statistical analysis
msec <- fun_summary(msec)
msec <- fun_analysis(msec,
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
fun_save(msec, "NCIH2286_LUNG", ".")

```

### Reproducibility

The source code is available in MicroSEC.R. 
Source data will be available at the Japanese 
Genotype-Phenotype Archive (http://trace.ddbj.nig.ac.jp/jga), 
which is hosted by the DNA Data Bank of Japan (the accession 
number is written in our paper). Each filtering process takes 
about 5–120 minutes per sample on a recommended machine, 
according to the depth and mutation amount of the sample.
