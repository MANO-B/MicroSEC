# MicroSEC: Microhomology-induced chimeric read-originating sequence error cleaning pipeline for FFPE samples
# 
# Version 1.0.3: parameter justification, version management with github (Oct. 29, 2020)
# Version 1.0.2: Indel length correction based on Bowtie2 penalty score (Oct. 28, 2020)
# Version 1.0.1: Correspond to the mouse genome. Homopolymer should have >=12 bases (Oct. 25, 2020)
# Version 1.0.0: Original version (Oct. 23, 2020)
# 
# Author: "Masachika Ikegami"
# Date: "Oct. 28, 2020"
# 
# This pipeline is designed for filtering mutations found in formalin-fixed and paraffin-embedded (FFPE) samples.
# The MicroSEC filter utilizes a statistical analysis, and the results for mutations with less than 10 supporting reads are not reliable.
# 
# Three files are nessesary for the analysis: mutation information file, BAM file, and mutation supporting read ID information file.
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
# This pipeline contains 8 filtering processes.
#
# Filter 1  : Shorter-supporting lengths distribute too short to occur (1-1 and 1-2).  
# Filter 1-1: P-values are less than the threshold_p(default: 10^(-6)).  
# Filter 1-2: The longest shorter-supporting lengths is shorter than 40% of the read length.  
# Filter 2  : Palindromic sequences exist within 150 bases (2-1 or 2-2).  
# Filter 2-1: A reverse complementary sequence consisting >= 30 bases exists on the opposite strand.  
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
# Rscript MicroSEC.R [working/output directory] [sample name] [mutation information excel file] [BAM file] [read ID information directory] [read length] [adaptor sequence] [progress bar Y/N] [sample type: Human or Mouse]
#
# Example
# Rscript MicroSEC.R /mnt/HDD8TB/FFPE SL_0002_L_FFPE_11-B /mnt/HDD8TB/FFPE/source/SL_mutation_somatic_200909.xlsx /mnt/HDD8TB/FFPE/source/SL_0002_L_FFPE_11_TDv4.realigned.bam /mnt/HDD8TB/FFPE/source/SL_mut_call_FFPE 150 AGATCGGAAGAGC Y Human
# Rscript MicroSEC.R /mnt/HDD8TB/FFPE SL_1010-N6-B /mnt/HDD8TB/FFPE/source/SL_mutation_somatic_200907.xlsx /mnt/HDD8TB/FFPE/source/SL_1010-N6_TDv4tw.realigned.bam /mnt/HDD8TB/FFPE/source/SL_1010-N6-B 124 AGATCGGAAGAGC Y Human
# Rscript MicroSEC.R /mnt/HDD8TB/FFPE PC9 /mnt/HDD8TB/FFPE/source/CCLE.xlsx /mnt/HDD8TB/FFPE/source/Cell_line/PC9_FEPE_Ag_TDv4.realigned.bam /mnt/HDD8TB/FFPE/source/PC9_FEPE 127 AGATCGGAAGAGC Y Human
#
#
#
# Sample name is set to the sample of interest in the Sample column of the mutation information file.
# Confirm the read length in the plathome
# Confirm the adaptor sequence; Todai Onco Panel ("AGATCGGAAGAGC")
# If you want to know the progress visually, [progress bar Y/N] should be Y.
# If you want to use the Filter 6, [GGGenome search Y/N] should be Y (requiring the internet access).
#
# Results are saved in a excel file.
# The explatation of the results is written in detail in the second sheet of the excel file.

# load necessary packages
library(tidyr)
library(openxlsx)
library(data.table)
library(R.utils)
library(stringr)
library(magrittr)
library(dplyr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicAlignments)
library(Rsamtools)
library(tcltk)
library(gtools)

# set arguments
args = commandArgs(trailingOnly = T)
wd = args[1]
SAMPLE_NAME = args[2]
MUTATION_FILE = args[3]
BAM_FILE = args[4]
MUTATION_SUPPORTING_READ_LIST = args[5]
READ_length = as.integer(args[6])
ADAPTOR_SEQ = args[7]
PROGRESS_BAR = args[8]
GENOME = args[9]

setwd(wd)

## Setting for follwing analyses

# hyperparameters
Width = 150
Pre_search_length_default = 0
Post_search_length_default = 20
Pre_search_length = Pre_search_length_default
Post_search_length = Post_search_length_default
Hairpin_search_length = 50
Minimum_Hairpin_length = 15
Short_Homology_search_length = 4
Minimum_Homology_search_length = 40
neighbor_length = 20
Homopolymer_length = 15
Laxness = 1
Length_Flag = 0
Max_mutation_search = 50

threshold_p = 10^(-6)
threshold_hairpin_ratio = 0.50
threshold_hairpin_length = 30
threshold_soft_clip_ratio = 0.90
threshold_short_length = 0.8
threshold_distant_homology = 0.2
threshold_low_quality_rate = 0.1

# initial value
CHROM = ""
MicroSEC = NULL
Homology_search = NULL
Pre_search_length = 0
Post_search_length = 0
Peri_seq_1 = ""
Peri_seq_2 = ""

# functions
fun_setting = function(PRE = Pre_search_length_default, POST = Post_search_length_default){
  Pre_search_length <<- PRE
  Post_search_length <<- POST
  Peri_seq_1 <<- substring(Neighbor_seq, neighbor_length + 1 - Pre_search_length, neighbor_length + Post_search_length + Alt_length)
  Peri_seq_2 <<- substring(Neighbor_seq, neighbor_length + 1 - Post_search_length, neighbor_length + Pre_search_length + Alt_length)
}

fun_hairpin_trimming = function(Hairpin_seq){
  if(mut_read_strand[[j]] == "+"){
    Hairpin_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ, subject = Hairpin_seq)
    Hairpin_seq = str_split(Hairpin_seq, ADAPTOR_SEQ)[[1]]
    if(length(Hairpin_seq)> 1){
      Hairpin_seq = Hairpin_seq[[1]]
    }
    Hairpin_seq = reverseComplement(DNAString(Hairpin_seq))
  }
  else{
    Hairpin_seq = reverseComplement(Hairpin_seq)
    Hairpin_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ, subject = Hairpin_seq)
    Hairpin_seq = str_split(Hairpin_seq, ADAPTOR_SEQ)[[1]]
    if(length(Hairpin_seq)> 1){
      Hairpin_seq = Hairpin_seq[[1]]
    }
    Hairpin_seq = DNAString(Hairpin_seq)
  }
  return(Hairpin_seq)
}

fun_hairpin_check = function(Hairpin_seq_tmp){
  Reverse_seq = matchPattern(pattern = Hairpin_seq_tmp, subject = Ref_seq, max.mismatch=0, min.mismatch=0,  with.indels=FALSE, fixed=TRUE)
  if(length(Reverse_seq) > 0){
    Hairpin_length <<- max(Hairpin_length, hair)
    return(1)
  }
  else{
    return(0)
  }
}

fun_repeat_check = function(Rep_A = DNAString(df_mutation[i,"Ref"]), Rep_B = DNAString(df_mutation[i,"Alt"]), Del = 1){
  Rep_B = Rep_B[2:nchar(Rep_B)]
  homo_tmp_1 = 0
  homo_tmp_2 = 0
  for(q in 1:nchar(Rep_B)){
    Rep_seq = Rep_B[q:nchar(Rep_B)]
    check_rep = TRUE
    Post_rep_status_tmp = 0
    for(r in (1 + Del * nchar(Rep_B)):(Width - 1)){
      if(check_rep &  Ref_seq[Width + r + 1] == Rep_seq[((r - 1) %% nchar(Rep_seq) + 1)]){
        Post_rep_status_tmp = Post_rep_status_tmp + 1
      }
      else{
        check_rep = FALSE
      }
    }
    if(Post_rep_status_tmp > 0){
      Post_rep_status <<- max(Post_rep_status, Post_rep_status_tmp)
      if(gsub(as.character(Rep_seq[1]), "", as.character(Rep_seq)) == ""){
        homo_tmp_1 = Post_rep_status_tmp + nchar(Rep_seq)
      }
    }
  }
  for(q in 1:nchar(Rep_B)){
    Rep_seq = Rep_B[1:q]
    check_rep = TRUE
    Pre_rep_status_tmp = 0
    for(r in 1:(Width + 1)){
      if(check_rep &  Ref_seq[Width + 2 - r] == Rep_seq[((nchar(Rep_seq) - r) %% nchar(Rep_seq) + 1)]){
        Pre_rep_status_tmp = Pre_rep_status_tmp + 1
      }
      else{
        check_rep = FALSE
      }
    }
    if(Pre_rep_status_tmp > 0){
      Pre_rep_status <<- max(Pre_rep_status, Pre_rep_status_tmp)
      if(gsub(as.character(Rep_seq[1]), "", as.character(Rep_seq)) == ""){
        homo_tmp_2 = Pre_rep_status_tmp + nchar(Rep_seq)
      }
    }
  }
  Homopolymer_status <<- max(homo_tmp_1, homo_tmp_2)
  if((Post_rep_status + Post_rep_status) > 0 & gsub(as.character(Rep_B[1]), "", as.character(Rep_B)) == ""){
    Homopolymer_status <<- homo_tmp_1 + homo_tmp_2 - nchar(Rep_B)
  }
}

zero = function(a,b)ifelse(a==0 & b==0, 0, a/b)

## Filtering

# load genomic sequence
if(GENOME == "Human"){
  genome = BSgenome.Hsapiens.UCSC.hg38
  Chr_No = 24
}
if(GENOME == "Mouse"){
  genome = BSgenome.Mmusculus.UCSC.mm10
  Chr_No = 21
}

# load somatic mutation list
df_mutation = read.xlsx(MUTATION_FILE, sheet = 1)
# load BAM file
file_BAM <- BamFile(BAM_FILE)
df_BAM = scanBam(file_BAM)[[1]]
df_BAM = df_BAM[c("qname", "rname", "strand", "cigar", "seq", "qual")]
# load mutation supporting read list
files = data.frame(file = list.files(MUTATION_SUPPORTING_READ_LIST, pattern=".gz", full.names=T))
df_mut_call = NULL
for(i in files$file){
  df_mut_call  = rbind(df_mut_call, read.csv(as.character(i), header=TRUE, sep="\t"))
}
# data formatting
if("Transition" %in% colnames(df_mutation)){
  df_mutation = df_mutation %>%
    dplyr::select(Sample, Gene, HGVS.p, Mut_type, `Total_QV>=20`, `%Alt`, Chr, Pos, Ref, Alt, SimpleRepeat_TRF, Neighborhood_sequence, Transition)
}
if(!"Transition" %in% colnames(df_mutation)){
  df_mutation = df_mutation %>%
    dplyr::select(Sample, Gene, HGVS.p, Mut_type, `Total_QV>=20`, `%Alt`, Chr, Pos, Ref, Alt, SimpleRepeat_TRF, Neighborhood_sequence)
  df_mutation$Transition = "NA"
}
df_mutation$Pos = as.integer(df_mutation$Pos)
df_mutation = df_mutation %>%
  dplyr::filter(Sample == SAMPLE_NAME) %>%
  dplyr::mutate(Ref = toupper(Ref)) %>%
  dplyr::mutate(Alt = toupper(Alt)) %>% 
  dplyr::mutate(Neighborhood_sequence = toupper(Neighborhood_sequence))
df_mut_call = df_mut_call %>%
  dplyr::select(Chr, Pos, Ref, Alt, Mut_ID, Mut) %>%
  dplyr::mutate(Ref = toupper(Ref)) %>%
  dplyr::mutate(Alt = toupper(Alt)) %>% 
  dplyr::mutate(Mut = toupper(Mut))
# analyze each somatic mutation
for(i in 1:length(df_mutation[,1])){
  Error = 0
  indel_status = 0
  indel_length = 0
  WITH_INDEL_1 = FALSE
  WITH_INDEL_2 = FALSE
  mut_call = logical(0)
  # extract mutation supporting reads
  if(df_mutation[i,"Chr"] != CHROM){
    CHROM = df_mutation[i,"Chr"]
    ID_No = df_BAM$rname == CHROM
    df_BAM_qname = df_BAM$qname[ID_No]
    df_BAM_seq = df_BAM$seq[ID_No]
    df_BAM_strand = df_BAM$strand[ID_No]
    df_BAM_cigar = df_BAM$cigar[ID_No]
    df_BAM_qual = df_BAM$qual[ID_No]
  }
  mut_read = df_mut_call %>% 
    dplyr::filter(Chr == df_mutation[i,"Chr"] & Pos == df_mutation[i,"Pos"])
  mut_type = str_split(df_mutation[i,"Mut_type"],"-")[[1]][[2]]
  if(mut_type == "ins"){
    indel_length = nchar(df_mutation[i,"Alt"]) - 1
  }
  if(mut_type == "del"){
    indel_length = nchar(df_mutation[i,"Ref"]) - 1
  }
  if(mut_type %in% c("ins", "del")){
    indel_status = 1
    WITH_INDEL_1 = TRUE
    WITH_INDEL_2 = TRUE
  }
  if(dim(mut_read)[[1]] == 0 & indel_status == 1){
    for(tmp in 1:Max_mutation_search){
      if(dim(mut_read)[[1]] == 0){
        mut_read = df_mut_call %>%
          dplyr::filter(Chr == df_mutation[i,"Chr"] & Pos == (df_mutation[i,"Pos"] - tmp))
        Error = tmp
      }
    }
  }
  if(dim(mut_read)[1] > 0){
    mut_detail = str_split(mut_read$Mut, pattern = ";")[[1]]
    mut_read_ID_list = str_split(mut_read$Mut_ID, pattern = ";")[[1]]
    if(indel_status == 0){
      mut_call = which(mut_detail == str_sub(df_mutation[i,"Alt"], start = 1, end = 1))
    }
    else if(mut_type == "del"){
      mut_call = which(mut_detail == str_replace(df_mutation[i,"Ref"], pattern=df_mutation[i,"Alt"], replacement=".-"))
    }
    else if(mut_type == "ins"){
      mut_call = which(mut_detail == str_replace(df_mutation[i,"Alt"], pattern=df_mutation[i,"Ref"], replacement=".+"))
    }
  }
  # if mutation supporting reads exist
  if(length(mut_call) > 0){
    Pre_support_length = 0
    Post_support_length = 0
    short_support_length = 0
    Pre_Minimum_length = Width
    Post_Minimum_length = Width
    Pre_rep_status = 0
    Post_rep_status = 0
    Homopolymer_status = 0
    FLAG_Hairpin = 0
    Hairpin_length = 0
    Soft_Clipped_read = 0
    Total_read = 0
    Low_quality_base_rate_under_Q18 = 0
    Homology_search_tmp = NULL
    distant_homology = 0
    mutated_ID = mut_read_ID_list[mut_call]
    mut_read_ID = str_sub(str_split(mutated_ID,",")[[1]], start=1, end=-2)
    mut_read_strand = str_sub(str_split(mutated_ID,",")[[1]], start=-1, end=-1)
    mut_read_strand = gsub("r", "-", gsub("f", "+", mut_read_strand))
    Ref_seq = genome[[df_mutation[i,"Chr"]]][(df_mutation[i,"Pos"] - Error - Width):(df_mutation[i,"Pos"] - Error + Width)]
    Neighbor_seq = df_mutation[i,"Neighborhood_sequence"]
    Alt_length = nchar(df_mutation[i,"Alt"])
    # sequence information around the mutation position
    fun_setting(PRE = Pre_search_length_default, POST = Post_search_length_default)
    near_list_1 = df_mut_call %>%
        dplyr::filter(Chr == df_mutation[i,"Chr"] &
                      Pos >= (df_mutation[i,"Pos"] - Error - Pre_search_length) &
                      Pos <= (df_mutation[i,"Pos"] + Post_search_length) &
                      Pos != (df_mutation[i,"Pos"] - Error))
    near_list_2 = df_mut_call %>%
        dplyr::filter(Chr == df_mutation[i,"Chr"] &
                      Pos >= (df_mutation[i,"Pos"] - Error - Post_search_length) &
                      Pos <= (df_mutation[i,"Pos"] + Pre_search_length) &
                      Pos != (df_mutation[i,"Pos"] - Error))
    Mut_near_1 = dim(near_list_1)[1]
    Mut_near_2 = dim(near_list_2)[1]
    if(all(nchar(near_list_1$Ref) != nchar(near_list_1$Alt))){
      WITH_INDEL_1 = TRUE
    }
    if(all(nchar(near_list_2$Ref) != nchar(near_list_2$Alt))){
      WITH_INDEL_2 = TRUE
    }
    if(PROGRESS_BAR == "Y"){
      cat(paste(" ", i, " / ", dim(df_mutation)[[1]], sep=""))
    }
    # short repeat around indel mutations
    if(indel_status == 1){
      if(mut_type == "ins"){
        fun_repeat_check(DNAString(df_mutation[i,"Ref"]), DNAString(df_mutation[i,"Alt"]), 0)
      }
      if(mut_type == "del"){
        fun_repeat_check(DNAString(df_mutation[i,"Alt"]), DNAString(df_mutation[i,"Ref"]), 1)
      }
    }
    # analyze each mutation supporting read
    if(PROGRESS_BAR == "Y"){
      pb = txtProgressBar(min = 0, max = max(1, length(mut_read_ID)), style = 3)
    }
    for(j in 1:length(mut_read_ID)){
      if(Length_Flag == 1){
        fun_setting(PRE = Pre_search_length_default, POST = Post_search_length_default)
      }
      Length_Flag = 0
      mut_position = 0
      FLAG_1 = Pre_search_length
      FLAG_2 = Pre_search_length
      # progress bar
      if(PROGRESS_BAR == "Y"){
        setTxtProgressBar(pb, j)
      }
      # specific read selection
      ID_No = df_BAM_qname == mut_read_ID[[j]]
      df_seq = df_BAM_seq[ID_No]
      df_strand = df_BAM_strand[ID_No]
      df_cigar = df_BAM_cigar[ID_No]
      df_qual = df_BAM_qual[ID_No]
      df_seq = df_seq[df_strand == mut_read_strand[[j]]]
      df_cigar = df_cigar[df_strand == mut_read_strand[[j]]]
      df_qual = df_qual[df_strand == mut_read_strand[[j]]]
      if(length(df_seq) > 1){
        df_cigar = df_cigar[which.max(width(df_seq))]
        df_qual = df_qual[which.max(width(df_seq)),]
        df_seq = df_seq[which.max(width(df_seq)),]
      }
      df_seq = df_seq[[1]]
      df_qual = as.vector(asc(as.character(df_qual[1])))
      if(length(df_seq) > 0){
        # determine mutation position in each read
        mutation_supporting_1 = matchPattern(Peri_seq_1, df_seq, max.mismatch=Mut_near_1 + Laxness, min.mismatch=0, with.indels=WITH_INDEL_1, fixed=FALSE)
        mutation_supporting_2 = matchPattern(Peri_seq_2, df_seq, max.mismatch=Mut_near_2 + Laxness, min.mismatch=0, with.indels=WITH_INDEL_2, fixed=FALSE)
        if(length(mutation_supporting_1) != 1 & length(mutation_supporting_2) != 1){
          Length_Flag = 1
          if(length(mutation_supporting_1) == 0 & length(mutation_supporting_2) == 0){
            for(Lax_1 in seq(3, 9, length=3)){
              for(Lax_2 in 2:0){
                if(length(mutation_supporting_1) != 1 & length(mutation_supporting_2) != 1){
                  fun_setting(PRE = Pre_search_length_default, POST = Post_search_length_default - Lax_1 * Laxness)
                  mutation_supporting_1 = matchPattern(Peri_seq_1, df_seq, max.mismatch=Mut_near_1 + Lax_2 * Laxness, min.mismatch=0, with.indels=WITH_INDEL_1, fixed=FALSE)
                  mutation_supporting_2 = matchPattern(Peri_seq_2, df_seq, max.mismatch=Mut_near_2 + Lax_2 * Laxness, min.mismatch=0, with.indels=WITH_INDEL_2, fixed=FALSE)
                }
              }
            }
          }
          if(length(mutation_supporting_1) != 1 & length(mutation_supporting_2) != 1){
            for(Lax_1 in seq(5, 15, length=3)){
              for(Lax_2 in 3:0){
                if(length(mutation_supporting_1) != 1 & length(mutation_supporting_2) != 1){
                  fun_setting(PRE = Pre_search_length_default + Lax_1 * Laxness, POST = Post_search_length_default)
                  mutation_supporting_1 = matchPattern(Peri_seq_1, df_seq, max.mismatch=Mut_near_1 + Lax_2 * Laxness, min.mismatch=0, with.indels=WITH_INDEL_1, fixed=FALSE)
                  mutation_supporting_2 = matchPattern(Peri_seq_2, df_seq, max.mismatch=Mut_near_2 + Lax_2 * Laxness, min.mismatch=0, with.indels=WITH_INDEL_2, fixed=FALSE)
                }
              }
            }
          }
        }
        if(length(mutation_supporting_1) == 1){
          mut_position = min(nchar(df_seq), start(mutation_supporting_1) + Pre_search_length)
          FLAG_1 = Hairpin_search_length
        }
        else if(length(mutation_supporting_2) == 1){
          mut_position = min(nchar(df_seq), start(mutation_supporting_2) + Post_search_length)
          FLAG_2 = Hairpin_search_length
        }
        if(mut_position > 0){
          Pre_support_length_tmp = 0
          Post_support_length_tmp = 0
          Total_read = Total_read + 1
          Reverse_seq = logical(0)
          FLAG_Hairpin_tmp = 0
          check_hairpin = 1
          # hairpin length calculation
          Hairpin_seq = fun_hairpin_trimming(df_seq[max(1,(mut_position - FLAG_1)):min(length(df_seq), (mut_position + FLAG_2 + Alt_length - 1))])
          if(Minimum_Hairpin_length < nchar(Hairpin_seq)){
            for(hair in Minimum_Hairpin_length:nchar(Hairpin_seq)){
              if(check_hairpin == 1){
                check_hairpin = fun_hairpin_check(Hairpin_seq_tmp = Hairpin_seq[1:hair])
                FLAG_Hairpin_tmp = max(FLAG_Hairpin_tmp, check_hairpin)
              }
            }
          }
          if(FLAG_Hairpin_tmp == 0){
            check_hairpin = 1
            Hairpin_seq = fun_hairpin_trimming(df_seq[max(1,(mut_position - Minimum_Hairpin_length)):min(length(df_seq), (mut_position + Minimum_Hairpin_length + Alt_length - 1))])
            if(Minimum_Hairpin_length < nchar(Hairpin_seq)){
              for(hair in Minimum_Hairpin_length:nchar(Hairpin_seq)){
                if(check_hairpin == 1){
                  check_hairpin = fun_hairpin_check(Hairpin_seq_tmp = Hairpin_seq[1:hair])
                  check_hairpin = max(check_hairpin, fun_hairpin_check(Hairpin_seq_tmp = Hairpin_seq[(nchar(Hairpin_seq) - hair + 1):nchar(Hairpin_seq)]))
                  FLAG_Hairpin_tmp = max(FLAG_Hairpin_tmp, check_hairpin)
                }
              }
            }
          }
          FLAG_Hairpin = FLAG_Hairpin + FLAG_Hairpin_tmp
          
          # supporting length calculation
          cigar_num = as.integer(str_split(df_cigar, "[:upper:]")[[1]])
          cigar_type = str_split(df_cigar, "[:digit:]+")[[1]][-1]
          tmp_pos = 1
          tmp_pos_indel = 0
          Soft_Clipped_read_tmp = 0
          for(k in 1:length(cigar_type)){
            if(cigar_type[k] == "H"){
              tmp_pos = tmp_pos + cigar_num[k]
            }
            if(cigar_type[k] == "I"){
              tmp_pos = tmp_pos + cigar_num[k]
            }
            if(cigar_type[k] == "S"){
              clipped_seq = df_seq[tmp_pos:(tmp_pos + cigar_num[k] - 1)]
              if(mut_read_strand[[j]] == "+"){
                clipped_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ, subject = clipped_seq)
                clipped_seq = str_split(clipped_seq, ADAPTOR_SEQ)[[1]]
                if(length(clipped_seq)> 1){
                  clipped_seq = clipped_seq[[1]]
                }
                clipped_seq = DNAString(clipped_seq)
              }
              else{
                clipped_seq = reverseComplement(clipped_seq)
                clipped_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ, subject = clipped_seq)
                clipped_seq = str_split(clipped_seq, ADAPTOR_SEQ)[[1]]
                if(length(clipped_seq)> 1){
                  clipped_seq = clipped_seq[[1]]
                }
                clipped_seq = DNAString(clipped_seq)
              }
              if(nchar(clipped_seq) > 0){
                Soft_Clipped_read_tmp = 1
              }
              tmp_pos = tmp_pos + cigar_num[k]
            }
            if(cigar_type[k] == "M"){
              if((tmp_pos + cigar_num[k] - 1) < mut_position){
                Pre_support_length_tmp = Pre_support_length_tmp + cigar_num[k]
              }
              else if(tmp_pos > (mut_position + Alt_length - 1)){
                Post_support_length_tmp = Post_support_length_tmp + cigar_num[k]
              }
              else{
                if(indel_status == 1){
                  Pre_support_length_tmp = Pre_support_length_tmp + (mut_position - tmp_pos)
                  Post_support_length_tmp = Post_support_length_tmp + (cigar_num[k] + tmp_pos - mut_position)
                }
                else{
                  Pre_support_length_tmp = Pre_support_length_tmp + (mut_position - tmp_pos)
                  Post_support_length_tmp = Post_support_length_tmp + (cigar_num[k] + tmp_pos - mut_position - Alt_length)
                }
              }
              tmp_pos = tmp_pos + cigar_num[k]
            }
          }
          
          # adjustment for consecutive snv
          if(mut_type == "snv"){
            Post_support_length_tmp = Post_support_length_tmp - Alt_length + 1
          }
          
          # read quality check
          Low_quality_base_rate_under_Q18 = Low_quality_base_rate_under_Q18 + sum(df_qual < 51) / nchar(df_seq)
          
          # query DNA sequence from GGGenome
          if(FLAG_Hairpin_tmp == 0 & FLAG_Hairpin == 0){
            if(indel_status == 1){
              Pre_Homology_search_seq = df_seq[1:min(length(df_seq), mut_position + Short_Homology_search_length + Post_rep_status + Alt_length - 1)]
              Post_Homology_search_seq = df_seq[max(1, mut_position - Short_Homology_search_length - Pre_rep_status):length(df_seq)]
            }
            if(indel_status == 0){
              Pre_Homology_search_seq = df_seq[1:min(length(df_seq), mut_position + Short_Homology_search_length + Alt_length)]
              Post_Homology_search_seq = df_seq[max(1, mut_position - Short_Homology_search_length):length(df_seq)]
            }
            if(mut_read_strand[[j]] == "-"){
              Pre_Homology_search_seq = fun_hairpin_trimming(Pre_Homology_search_seq)
              Pre_Homology_search_seq = reverseComplement(Pre_Homology_search_seq)
            }
            if(mut_read_strand[[j]] == "+"){
              Post_Homology_search_seq = fun_hairpin_trimming(Post_Homology_search_seq)
              Post_Homology_search_seq = reverseComplement(Post_Homology_search_seq)
            }
            Homology_search_tmp = rbind(Homology_search_tmp, data.frame(SAMPLE_NAME, Chr = df_mutation[i,"Chr"],  Pos = df_mutation[i,"Pos"], Ref = df_mutation[i,"Ref"], Alt = df_mutation[i,"Alt"], Direction = "Pre", Seq = as.character(Pre_Homology_search_seq)))
            Homology_search_tmp = rbind(Homology_search_tmp, data.frame(SAMPLE_NAME, Chr = df_mutation[i,"Chr"],  Pos = df_mutation[i,"Pos"], Ref = df_mutation[i,"Ref"], Alt = df_mutation[i,"Alt"], Direction = "Post", Seq = as.character(Post_Homology_search_seq)))
          }

          # summary
          Soft_Clipped_read = Soft_Clipped_read + Soft_Clipped_read_tmp
          Pre_support_length = max(Pre_support_length, Pre_support_length_tmp)
          Post_support_length = max(Post_support_length, Post_support_length_tmp)
          short_support_length = max(short_support_length, min(Pre_support_length_tmp, Post_support_length_tmp))
          Pre_Minimum_length = min(Pre_Minimum_length, Pre_support_length_tmp)
          Post_Minimum_length = min(Post_Minimum_length, Post_support_length_tmp)
        }
      }
    }
    # data formatting
    MicroSEC_tmp = df_mutation[i,] %>% dplyr::mutate(
      READ_length = READ_length, 
      mut_type = mut_type,
      Alt_length = Alt_length, 
      Total_read = Total_read, 
      Soft_Clipped_read = Soft_Clipped_read, 
      FLAG_Hairpin = FLAG_Hairpin, 
      Hairpin_length = Hairpin_length, 
      Pre_support_length = Pre_support_length, 
      Post_support_length = Post_support_length, 
      short_support_length = short_support_length, 
      Pre_Minimum_length = Pre_Minimum_length, 
      Post_Minimum_length = Post_Minimum_length, 
      Low_quality_base_rate_under_Q18 = Low_quality_base_rate_under_Q18 / Total_read,
      Pre_rep_status = Pre_rep_status, 
      Post_rep_status = Post_rep_status, 
      Homopolymer_status = Homopolymer_status,
      indel_status = indel_status,
      indel_length = indel_length,
      distant_homology = distant_homology
    )
    if(FLAG_Hairpin == 0){
      Homology_search = rbind(Homology_search, Homology_search_tmp)
    }
  }
  else{
    MicroSEC_tmp = df_mutation[i,] %>% dplyr::mutate(
      READ_length = READ_length, 
      mut_type = mut_type,
      Alt_length = 0, 
      Total_read = 0, 
      Soft_Clipped_read = 0, 
      FLAG_Hairpin = 0, 
      Hairpin_length = 0, 
      Pre_support_length = 0, 
      Post_support_length = 0, 
      short_support_length = 0, 
      Pre_Minimum_length = 0, 
      Post_Minimum_length = 0, 
      Low_quality_base_rate_under_Q18 = 0,
      Pre_rep_status = 0, 
      Post_rep_status = 0, 
      Homopolymer_status = 0,
      indel_status = 0,
      indel_length = 0,
      distant_homology = 0
    )
  }
  MicroSEC = rbind(MicroSEC, MicroSEC_tmp)
}

# search homologous sequences from the entire genome by matchPDict function
MicroSEC$distant_homology = 0
df_distant = Homology_search
df_distant = df_distant %>% dplyr::mutate(
  Seq = as.character(Seq)
)
df_distant = df_distant %>% dplyr::mutate(
  length = nchar(Seq)
)
df_distant$Number = floor(seq(1,(dim(df_distant)[1]/2 + 0.5),0.5))
MAX_No = max(df_distant$Number)
df_distant = df_distant[df_distant$length >= Minimum_Homology_search_length,]
df_distant = df_distant %>% dplyr::mutate(
  fixed_seq_Pre = str_sub(Seq, length - Minimum_Homology_search_length + 1, length),
  fixed_seq_Post = str_sub(Seq, 1, Minimum_Homology_search_length)
)
df_distant = df_distant %>% dplyr::mutate(
  fixed_seq = ifelse(Direction == "Pre", fixed_seq_Pre, fixed_seq_Post)
)
df_distant = df_distant[!str_detect(df_distant$fixed_seq, pattern="N"),]
Search_Seq_f = DNAStringSet(df_distant$fixed_seq)
Search_Seq_r = reverseComplement(Search_Seq_f)
Search_Seq_f = PDict(Search_Seq_f)
Search_Seq_r = PDict(Search_Seq_r)
distant_homology = rep(0, length(Search_Seq_f))
if(PROGRESS_BAR == "Y"){
  pb = txtProgressBar(min = 0, max = Chr_No, style = 3)
}
for(seqname in 1:Chr_No){
  if(PROGRESS_BAR == "Y"){
    setTxtProgressBar(pb, seqname)
  }
  target = genome[[seqname]]
  distant_homology = distant_homology + countPDict(Search_Seq_f, target)
  distant_homology = distant_homology + countPDict(Search_Seq_r, target)
}
df_distant$distant_homology = distant_homology
if(PROGRESS_BAR == "Y"){
  pb = txtProgressBar(min = 0, max = max(1, dim(df_distant)[1]), style = 3)
}
for(i in 1:MAX_No){
  if(PROGRESS_BAR == "Y"){
    setTxtProgressBar(pb, i)
  }
  if(dim(df_distant[df_distant$Number == i,])[1] == 1){
    tmp_distant = df_distant[df_distant$Number == i,][1,]
    if(tmp_distant$distant_homology > 0){
      MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME & MicroSEC$Chr == tmp_distant$Chr & MicroSEC$Pos == tmp_distant$Pos & MicroSEC$Ref == tmp_distant$Ref & MicroSEC$Alt == tmp_distant$Alt,]$distant_homology =  MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME & MicroSEC$Chr == tmp_distant$Chr & MicroSEC$Pos == tmp_distant$Pos & MicroSEC$Ref == tmp_distant$Ref & MicroSEC$Alt == tmp_distant$Alt,]$distant_homology + 1
    }
  }
  if(dim(df_distant[df_distant$Number == i,])[1] == 2){
    tmp_distant = df_distant[df_distant$Number == i,][1,]
    tmp_distant_2 = df_distant[df_distant$Number == i,][2,]
    if(tmp_distant$distant_homology > 0 | tmp_distant_2$distant_homology > 0){
      MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME & MicroSEC$Chr == tmp_distant$Chr & MicroSEC$Pos == tmp_distant$Pos & MicroSEC$Ref == tmp_distant$Ref & MicroSEC$Alt == tmp_distant$Alt,]$distant_homology =  MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME & MicroSEC$Chr == tmp_distant$Chr & MicroSEC$Pos == tmp_distant$Pos & MicroSEC$Ref == tmp_distant$Ref & MicroSEC$Alt == tmp_distant$Alt,]$distant_homology + 1
    }
  }
}


# statistical analysis
MicroSEC = MicroSEC %>% dplyr::mutate(
  distant_homology_rate = zero(distant_homology, Total_read)
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  Pre_Minimum_length = ifelse(indel_status == 1, (((Pre_Minimum_length - (indel_length + 1)) + abs(Pre_Minimum_length - (indel_length + 1))) / 2) + (indel_length + 1), Pre_Minimum_length),
  Post_Minimum_length = ifelse(indel_status == 1, (((Post_Minimum_length - (indel_length + 1)) + abs(Post_Minimum_length - (indel_length + 1))) / 2) + (indel_length + 1), Post_Minimum_length),
  Pre_support_length = ifelse(indel_status == 1, (((Pre_support_length - (READ_length - indel_length - 1)) - abs(Pre_support_length - (READ_length - indel_length - 1))) / 2) + (READ_length - indel_length - 1), Pre_support_length),
  Post_support_length = ifelse(indel_status == 1, (((Post_support_length - (READ_length - indel_length - 1)) - abs(Post_support_length - (READ_length - indel_length - 1))) / 2) + (READ_length - indel_length - 1), Post_support_length)
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  Pre_Minimum_length = ifelse(indel_status == 1, (((Pre_Minimum_length - Pre_rep_status) + abs(Pre_Minimum_length - Pre_rep_status)) / 2) + Pre_rep_status, Pre_Minimum_length),
  Post_Minimum_length = ifelse(indel_status == 1, (((Post_Minimum_length - Post_rep_status) + abs(Post_Minimum_length -Post_rep_status)) / 2) + Post_rep_status, Post_Minimum_length),
  Pre_support_length = ifelse(indel_status == 1, (((Pre_support_length - (READ_length - Post_rep_status)) - abs(Pre_support_length - (READ_length - Post_rep_status))) / 2) + (READ_length - Post_rep_status), Pre_support_length),
  Post_support_length = ifelse(indel_status == 1, (((Post_support_length - (READ_length - Pre_rep_status)) - abs(Post_support_length - (READ_length - Pre_rep_status))) / 2) + (READ_length - Pre_rep_status), Post_support_length)
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  shortest_support_length = (((Pre_Minimum_length - Post_Minimum_length) - abs(Pre_Minimum_length - Post_Minimum_length)) / 2) + Post_Minimum_length,
  minimum_length_1 = (((Pre_rep_status - (indel_length + 1)) + abs(Pre_rep_status - (indel_length + 1))) / 2) + (indel_length + 1),
  minimum_length_2 = (((Post_rep_status - (indel_length + 1)) + abs(Post_rep_status - (indel_length + 1))) / 2) + (indel_length + 1)
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  minimum_length_1 = ifelse(indel_status == 1, minimum_length_1, 0),
  minimum_length_2 = ifelse(indel_status == 1, minimum_length_2, 0)
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  minimum_length = (((minimum_length_1 -  minimum_length_2) - abs(minimum_length_1 -  minimum_length_2)) / 2) +  minimum_length_2
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  short_support_length = (((short_support_length - minimum_length) + abs(short_support_length - minimum_length)) / 2) + minimum_length,
  Altered_length = ifelse(mut_type == "snv", Alt_length, ifelse(mut_type == "ins", indel_length, 0)),
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  short_support_length_adjust = short_support_length - shortest_support_length + 1,
  Pre_support_length_adjust = Pre_support_length - Pre_Minimum_length + 1,
  Post_support_length_adjust = Post_support_length - Post_Minimum_length + 1,
  Half_length_adjust = as.integer((READ_length - Altered_length) / 2)  - minimum_length + 1,
  Total_length_adjust = READ_length - Altered_length - minimum_length_1 - minimum_length_2 + 1
)
MicroSEC = MicroSEC %>% dplyr::select(-minimum_length_1, -minimum_length_2, -minimum_length, -Altered_length, -shortest_support_length)
MicroSEC = MicroSEC %>% dplyr::mutate(
  Short_short_support = (short_support_length_adjust <= threshold_short_length * Half_length_adjust),
  Short_pre_support = (Pre_support_length_adjust <= threshold_short_length * Total_length_adjust),
  Short_post_support = (Post_support_length_adjust <= threshold_short_length * Total_length_adjust),
  High_rate_Q18 = ifelse((Low_quality_base_rate_under_Q18 < threshold_low_quality_rate), TRUE, FALSE),
  prob_Filter_1 = (short_support_length_adjust / Half_length_adjust) ^ Total_read,
  prob_Filter_3_pre = (Pre_support_length_adjust / Total_length_adjust) ^ Total_read,
  prob_Filter_3_post = (Post_support_length_adjust / Total_length_adjust) ^ Total_read
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  prob_Filter_1 = ifelse((prob_Filter_1 > 1), 1, prob_Filter_1),
  prob_Filter_3_pre = ifelse((prob_Filter_3_pre > 1), 1, prob_Filter_3_pre),
  prob_Filter_3_post = ifelse((prob_Filter_3_post > 1), 1, prob_Filter_3_post)
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  Filter_1_mutation_intra_hairpin_loop = ifelse((Short_short_support & prob_Filter_1 < threshold_p), TRUE, FALSE), 
  Filter_2_hairpin_structure = ifelse((zero(FLAG_Hairpin, Total_read) > threshold_hairpin_ratio | Hairpin_length >= threshold_hairpin_length), TRUE, FALSE),
  Filter_3_microhomology_induced_mutation = ifelse((High_rate_Q18 & ((prob_Filter_3_pre <= threshold_p & Short_pre_support) | (prob_Filter_3_post <= threshold_p & Short_post_support))), TRUE, FALSE),
  Filter_4_soft_clipping = ifelse((zero(Soft_Clipped_read, Total_read)) > threshold_soft_clip_ratio, TRUE, FALSE),
  Filter_5_highly_homologous_region = ifelse((distant_homology_rate >= threshold_distant_homology), TRUE, FALSE),
  Filter_6_simple_repeat = ifelse((SimpleRepeat_TRF == "Y"), TRUE, FALSE), 
  Filter_7_C_to_T_artifact = ifelse((Transition == "C>T_g_FFPE"), TRUE, FALSE),
  Filter_8_mutation_at_homopolymer = ifelse((Homopolymer_status  >= Homopolymer_length), TRUE, FALSE)
)
MicroSEC = MicroSEC %>% dplyr::mutate(
  MicroSEC_filter_1234 = Filter_1_mutation_intra_hairpin_loop | Filter_2_hairpin_structure | Filter_3_microhomology_induced_mutation | Filter_4_soft_clipping,
  MicroSEC_filter_12345 = Filter_1_mutation_intra_hairpin_loop | Filter_2_hairpin_structure | Filter_3_microhomology_induced_mutation | Filter_4_soft_clipping | Filter_5_highly_homologous_region,
  MicroSEC_filter_all = Filter_1_mutation_intra_hairpin_loop | Filter_2_hairpin_structure | Filter_3_microhomology_induced_mutation | Filter_4_soft_clipping | Filter_5_highly_homologous_region | Filter_6_simple_repeat | Filter_7_C_to_T_artifact | Filter_8_mutation_at_homopolymer
)

# explanation
MicroSEC_explain = data.frame(
  Alt_length  = "The length of altered sequence",
  Total_read = "The number of reads supporting the mutation: [Total_read >= 10] filtering is strongly recommended",
  mut_type = "Type of the mutation",
  Soft_Clipped_read = "The number of soft-clipped reads",
  FLAG_Hairpin = "The number of reads containing a reverse complementary sequence of the opposite strand consisting >= 15 bases",
  Hairpin_length = "The longest length of reverse complementary sequences existing on the opposite strand",
  Pre_support_length = "The longest mapped length of reads of the 5' sequence of the mutation",
  Post_support_length = "The longest mapped length of reads of the 3' sequence of the mutation",
  short_support_length = "The longest length of the shorter mapped length around the mutation of all mutation-supporting reads",
  Pre_Minimum_length = "The shortest mapped length of reads of the 5' sequence of the mutation",
  Post_Minimum_length = "The shortest mapped length of reads of the 3' sequence of the mutation",
  Low_quality_base_rate_under_Q18 = "mean rate of base quality under Q18",
  Pre_rep_status = "Adjustment for 5' Supporting lengths considering small repeat sequences around the mutations",
  Post_rep_status = "Adjustment for 3' Supporting lengths considering small repeat sequences around the mutations",
  Homopolymer_status = "The length of the homopolymer which the mutation located inside",
  indel_status = "Whether the type of the mutation is insertion/deletion or not",
  indel_length = "The length of insertion/deletion",
  distant_homology = "The number of chimeric reads comprising two distant regions",
  distant_homology_rate = "The rate of chimeric reads comprising two distant regions",
  short_support_length_adjust = "Adjusted short_support_length",
  Pre_support_length_adjust = "Adjusted Pre_support_length",
  Post_support_length_adjust = "Adjusted Post_support_length",
  Half_length_adjust = "Adjusted length of the half of the read length",
  Total_length_adjust = "Adjusted  length of the read length",
  Short_short_support = "short_support_length_adjust <= 80% of Half_length_adjust",
  Short_pre_support = "Pre_support_length_adjust <= 80% of Total_length_adjust",
  Short_post_support = "Post_support_length_adjust <= 80% of Total_length_adjust",
  High_rate_Q18 = "The rate of base quality under Q18 is <10% around the mutation" ,
  prob_Filter_1 = "The probability of the occurence of the shorter-supporting lengths distribution (threshold_p: 10^(-6))",
  prob_Filter_3_pre = "The probability of the occurence of the pre-supporting lengths distribution (threshold_p: 10^(-6))",
  prob_Filter_3_post = "The probability of the occurence of the post-supporting lengths distribution (threshold_p: 10^(-6))",
  Filter_1_mutation_intra_hairpin_loop = "Shorter-supporting lengths are too short to occur",
  Filter_2_hairpin_structure = "Palindromic sequences exist within 150 bases",
  Filter_3_microhomology_induced_mutation = "3’-/5’-supporting lengths are too densely distributed to occur",
  Filter_4_soft_clipping = ">=90% mutation-supporting reads are soft-clipped",
  Filter_5_highly_homologous_region = ">=20% mutations were called by chimeric reads comprising two distant regions",
  Filter_6_simple_repeat = "Mutations locating at simple repeat sequences",
  Filter_7_C_to_T_artifact = "C>T_g false positive calls in FFPE samples",
  Filter_8_mutation_at_homopolymer =  "Indel mutations located inside a >=15 homopolymer",
  MicroSEC_filter_1234 = "One or more filters of Fitler 1, 2, 3, or 4 are TRUE",
  MicroSEC_filter_12345 = "One or more filters of Fitler 1, 2, 3, 4, or 5 are TRUE",
  MicroSEC_filter_all= "One or more filters are TRUE"
)
MicroSEC_explain = data.frame(t(MicroSEC_explain))
colnames(MicroSEC_explain) = c("Explanation")
MicroSEC_explain_Name = data.frame(rownames(MicroSEC_explain))
colnames(MicroSEC_explain_Name) = c("Name")
MicroSEC_explain = cbind(MicroSEC_explain_Name, MicroSEC_explain)
colnames(MicroSEC_explain) = c("Name", "Explanation")

# save the results
NewWb = createWorkbook()
addWorksheet(wb = NewWb, sheetName = "MicroSEC_results", gridLines = TRUE)
addWorksheet(wb = NewWb, sheetName = "MicroSEC_explanation", gridLines = TRUE)
writeData(wb = NewWb, sheet = "MicroSEC_results", x = MicroSEC, xy = c(1,1), borders = "all", withFilter=TRUE)
writeData(wb = NewWb, sheet = "MicroSEC_explanation", x = MicroSEC_explain, xy = c(1,1), borders = "all")
freezePane(wb = NewWb, sheet = "MicroSEC_results", firstActiveRow = 2)
setColWidths(wb = NewWb, sheet = "MicroSEC_results", cols = 1:ncol(MicroSEC), widths = "auto")
setColWidths(wb = NewWb, sheet = "MicroSEC_explanation", cols = 1:ncol(MicroSEC_explain), widths = "auto")

saveWorkbook(wb = NewWb, file = paste(wd, "/MicroSEC_", SAMPLE_NAME, "_", Sys.Date(), ".xlsx", sep=""), overwrite = TRUE)
