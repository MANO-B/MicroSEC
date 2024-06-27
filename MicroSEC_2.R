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

Tumor = "BREAST"
wd = paste0("/tmp/MicroSEC/analysis_result/", Tumor)
dd = paste0("/tmp/MicroSEC/analysis_data/", Tumor)
sample_list_filename = paste0(wd, "/Sample_list.tsv")
analysis_cases = '/tmp/MicroSEC/clinical_data/all_cases.csv'
simple_repeat_list = "/tmp/MicroSEC/reference/simpleRepeat.bed.gz"
file_info_dir = paste0('/tmp/MicroSEC/clinical_data/', Tumor)
chromosomes = paste0("chr", c(seq_len(22),"X", "Y"))
progress_bar = "Y"
reference_genome = "/tmp/MicroSEC/reference/Homo_sapiens_assembly38.fasta"

options(digits.secs=5)
organism = "hg38"
fun_load_genome(organism)
fun_load_chr_no(organism)
setwd(wd)

print("Sample information loading...")

if(!dir.exists("/tmp/MicroSEC/analysis_result")){
  dir.create("/tmp/MicroSEC/analysis_result")
}
if(!dir.exists(wd)){
  dir.create(wd)
}
if(!dir.exists(paste0(wd, "/MicroSEC_result"))){
  dir.create(paste0(wd, "/MicroSEC_result"))
}

if(file.exists(sample_list_filename)){
  sample_info = read.csv(sample_list_filename, header=FALSE, stringsAsFactors = FALSE, sep = "\t")
} else{
  sample_list = read.csv(analysis_cases, header=TRUE, stringsAsFactors = FALSE, sep = ",")
  
  sample_list_F1_tumor = sample_list %>%
    dplyr::filter(scanCantypcat == Tumor & biospec.bescatEn == "FoundationOne CDx") %>%
    dplyr::select(hashregistrationid)
  colnames(sample_list_F1_tumor) = "sample_name"
  sample_list_F1_tumor$mutation_file = paste0(dd,"/",sample_list_F1_tumor$sample_name,"/mutation.tsv.gz")
  sample_list_F1_tumor$bam_file = paste0(dd,"/",sample_list_F1_tumor$sample_name,"/cram/",str_replace_all(sample_list_F1_tumor$sample_name, "-", ""),".markdup.cram")
  sample_list_F1_tumor$read_length = 49
  sample_list_F1_tumor$adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  sample_list_F1_tumor$adapter2 = "ACATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  sample_list_F1_tumor$refgenome = "hg38"
  sample_list_F1_tumor$panel = "F1"
  
  sample_list_NOP_tumor = sample_list %>%
    dplyr::filter(scanCantypcat == Tumor & biospec.bescatEn == "NCC OncoPanel") %>%
    dplyr::select(hashregistrationid)
  colnames(sample_list_NOP_tumor) = "sample_name"
  sample_list_NOP_tumor$mutation_file = paste0(dd,"/",sample_list_NOP_tumor$sample_name,"/mutation_tumor.tsv.gz")
  sample_list_NOP_tumor$bam_file = paste0(dd,"/",sample_list_NOP_tumor$sample_name,"/cram/",str_replace_all(sample_list_NOP_tumor$sample_name, "-", ""),"_tumor.markdup.cram")
  sample_list_NOP_tumor$read_length = 151
  sample_list_NOP_tumor$adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  sample_list_NOP_tumor$adapter2 = "ACATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  sample_list_NOP_tumor$refgenome = "hg38"
  sample_list_NOP_tumor$panel = "NOP_T"
  
  sample_list_NOP_normal = sample_list %>%
    dplyr::filter(scanCantypcat == Tumor & biospec.bescatEn == "NCC OncoPanel") %>%
    dplyr::select(hashregistrationid)
  colnames(sample_list_NOP_normal) = "sample_name"
  sample_list_NOP_normal$mutation_file = paste0(dd,"/",sample_list_NOP_normal$sample_name,"/mutation_normal.tsv.gz")
  sample_list_NOP_normal$bam_file = paste0(dd,"/",sample_list_NOP_normal$sample_name,"/cram/",str_replace_all(sample_list_NOP_normal$sample_name, "-", ""),"_normal.markdup.cram")
  sample_list_NOP_normal$read_length = 151
  sample_list_NOP_normal$adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  sample_list_NOP_normal$adapter2 = "ACATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  sample_list_NOP_normal$refgenome = "hg38"
  sample_list_NOP_normal$panel = "NOP_N"
  
  sample_info = rbind(sample_list_F1_tumor, sample_list_NOP_tumor, sample_list_NOP_normal)
  rm(sample_list, sample_list_F1_tumor, sample_list_NOP_tumor, sample_list_NOP_normal)
  sample_info = sample_info %>% dplyr::arrange(sample_name)
  write_tsv(progress = F, sample_info,file = sample_info_filename,col_names = FALSE)
}

names = c("sample_name","mutation_file","bam_file","read_length","adapter1","adapter2","refgenome", "panel")
colnames(sample_info) = names

for(sample in seq_len(dim(sample_info)[1])){
  sample_name = sample_info[sample, 1]
  mutation_file = sample_info[sample, 2]
  bam_file = sample_info[sample, 3]
  read_length = sample_info[sample, 4]
  adapter_1 = sample_info[sample, 5]
  adapter_2 = sample_info[sample, 6]
  organism = sample_info[sample, 7]
  panel = sample_info[sample, 8]
  
  start_timing = format(now(), "%H:%M:%OS")
  if(!file.exists(sample_name)){
    write(start_timing, append = F, sample_name)
    Sys.sleep(0.01)
    if(scan(sample_name, what=character()) == start_timing){
      write("Analyzing now...", append = F, sample_name)
      print(paste0("Sample No: ", sample, "/", dim(sample_info)[1]))
      if(panel %in% c("F1", "NOP_T", "NOP_N")){
        file_url_info = read_csv(col_names = FALSE, progress = F, show_col_types = F,paste0(file_info_dir, "/", Tumor, "_files_", floor(sample/10), ".csv"))
        file_crai = paste0(bam_file,".crai")
        if(panel == "F1"){
          url_cram = (file_url_info %>% dplyr::filter(X1 == sample_name & X5 == "cram" & !str_detect(X7, "crai")))$X9
        }
        if(panel == "NOP_T"){
          url_cram = (file_url_info %>% dplyr::filter(X1 == sample_name & X5 == "cram" &
                                                        str_detect(X7, "tumor") & !str_detect(X7, "crai")))$X9
        }
        if(panel == "NOP_N"){
          url_cram = (file_url_info %>% dplyr::filter(X1 == sample_name & X5 == "cram" &
                                                        str_detect(X7, "normal") & !str_detect(X7, "crai")))$X9
        }
        rm(file_url_info)
        list_exist = FALSE
        output = paste0(wd, "/MicroSEC_result/", sample_name,"_", panel, ".tsv.gz")
        bam_file_tmp = paste0(bam_file, ".tmp.bam")
        bam_file_slim = paste0(bam_file, ".SLIM.bam")
        tmp = data.frame(0)
        print("Preparing mutation information file...")
        if(!file.exists(output)){
          if(panel == "NOP_N"){
            mutation_raw_file = paste0(dd,"/",sample_name,
                                       "/germline_short_variant/",
                                       str_replace_all(sample_name, "-", ""),
                                       "_normal_rare-variant.vep.txt")
            tmp = read_tsv(progress = F, show_col_types = F,mutation_raw_file, skip = 1)
            if(dim(tmp)[1]!=0){
              mutation_table = data.frame(tmp)
              df_mutation = mutation_table %>% dplyr::select(
                Consequence, EXON, INTRON, Protein_position, Amino_acids, CLIN_SIG,
                gnomADg, gnomADg_AF, ToMMo, ToMMo_AF, GENIE_Mut, GENIE_Mut_sample_count)
              df_mutation$Sample = sample_name
              df_mutation$Chr = mutation_table$CHROM
              df_mutation$Pos = mutation_table$POS
              df_mutation$Ref = mutation_table$REF
              df_mutation$Alt = mutation_table$ALT
              df_mutation$Mismatch_Rate_Tumor = floor(1000*mutation_table$AD_FORMAT_ALT / mutation_table$DP_FORMAT)/1000
              df_mutation$Strand_Bias_Tumor = 0
              df_mutation$EBCall_Pval = 0
              df_mutation$HGVS.c = mutation_table$HGVSc
              df_mutation$HGVS.p = mutation_table$HGVSp
              df_mutation$Gene = mutation_table$SYMBOL
              df_mutation$Mut_type = "-"
              df_mutation$`Total_QV>=20` = 0
              df_mutation$`%Alt` = 0
              df_mutation$SimpleRepeat_TRF = "N"
              df_mutation$Neighborhood_sequence = "-"
              df_mutation$Transition = "-"
              df_mutation$Panel = "NOP_N"
              
              
              simple_repeat_info = data.frame(read_tsv(progress = F, col_names = FALSE, show_col_types = F,simple_repeat_list))[,1:3]
              simple_repeat_info = simple_repeat_info %>% filter(X1 %in% chromosomes)
              chr_now = "chr1"
              simple_repeat_now = simple_repeat_info %>% filter(X1 == chr_now)
              for(k in seq_len(length(df_mutation$SimpleRepeat_TRF))){
                if(df_mutation$Chr[k] != chr_now){
                  chr_now = df_mutation$Chr[k]
                  simple_repeat_now = simple_repeat_info %>% filter(X1 == chr_now)
                }
                if(length((simple_repeat_now %>% dplyr::filter(X1 == df_mutation$Chr[k] &
                                                               X2 <= df_mutation$Pos[k] &
                                                               X3 >= df_mutation$Pos[k]))$X1)){
                  df_mutation$SimpleRepeat_TRF[k] = "Y"
                }
              }
              rm(simple_repeat_info, simple_repeat_now, mutation_table)
            }
          }
          if(panel == "NOP_T"){
            mutation_raw_file = paste0(dd,"/",sample_name,
                                       "/somatic_short_variant/",
                                       str_replace_all(sample_name, "-", ""),
                                       "_tumor_gcatfilter_mutation.filt.vep.txt")
            tmp = read_tsv(progress = F, show_col_types = F,mutation_raw_file, skip = 1)
            if(dim(tmp)[1]!=0){
              mutation_table = data.frame(tmp)
              df_mutation = mutation_table %>% dplyr::select(
                Chr, Pos, Ref, Alt, Mismatch_Rate_Tumor, Strand_Bias_Tumor, EBCall_Pval,
                Consequence, EXON, INTRON, Protein_position, Amino_acids, CLIN_SIG,
                gnomADg, gnomADg_AF, ToMMo, ToMMo_AF, GENIE_Mut, GENIE_Mut_sample_count)
              df_mutation$Sample = sample_name
              df_mutation$HGVS.c = mutation_table$HGVSc
              df_mutation$HGVS.p = mutation_table$HGVSp
              df_mutation$Gene = mutation_table$SYMBOL
              df_mutation$Mut_type = "-"
              df_mutation$`Total_QV>=20` = 0
              df_mutation$`%Alt` = 0
              df_mutation$SimpleRepeat_TRF = mutation_table$Simple_Repeat_Region
              df_mutation$Neighborhood_sequence = "-"
              df_mutation$Transition = "-"
              df_mutation$Panel = "NOP_T"
              rm(mutation_table)
            }
          }
          if(panel == "F1"){
            mutation_raw_file = paste0(dd,"/",sample_name,
                                       "/somatic_short_variant/",
                                       str_replace_all(sample_name, "-", ""),
                                       "_gcatfilter_mutation.filt.vep.txt")
            tmp = read_tsv(progress = F, show_col_types = F,mutation_raw_file, skip = 1)
            if(dim(tmp)[1]!=0){
              mutation_table = data.frame(tmp)
              df_mutation = mutation_table %>% dplyr::select(
                Chr, Pos, Ref, Alt, Mismatch_Rate_Tumor, Strand_Bias_Tumor, EBCall_Pval,
                Consequence, EXON, INTRON, Protein_position, Amino_acids, CLIN_SIG,
                gnomADg, gnomADg_AF, ToMMo, ToMMo_AF, GENIE_Mut, GENIE_Mut_sample_count)
              df_mutation$Sample = sample_name
              df_mutation$HGVS.c = mutation_table$HGVSc
              df_mutation$HGVS.p = mutation_table$HGVSp
              df_mutation$Gene = mutation_table$SYMBOL
              df_mutation$Mut_type = "-"
              df_mutation$`Total_QV>=20` = 0
              df_mutation$`%Alt` = 0
              df_mutation$SimpleRepeat_TRF = mutation_table$Simple_Repeat_Region
              df_mutation$Neighborhood_sequence = "-"
              df_mutation$Transition = "-"
              df_mutation$Panel = "F1"
              rm(mutation_table)
            }
          }
          if(dim(tmp)[1]!=0){
            df_mutation = df_mutation %>% dplyr::mutate(
              SimpleRepeat_TRF = case_when(
                is.na(SimpleRepeat_TRF) ~ "N",
                SimpleRepeat_TRF == "Y" ~ "Y",
                SimpleRepeat_TRF == "N" ~ "N",
                SimpleRepeat_TRF == TRUE ~ "Y",
                TRUE ~ "Y"
              )
            )
            df_mutation = df_mutation %>% dplyr::mutate(
              Alt = case_when(
                Alt == "*" ~ str_sub(Ref, 1, 1),
                TRUE ~ Alt
              )
            )
            df_mutation$RefLen = nchar(df_mutation$Ref)
            df_mutation$AltLen = nchar(df_mutation$Alt)
            df_mutation = df_mutation %>% dplyr::mutate(
              Mut_type = case_when(
                RefLen == AltLen ~ paste0(RefLen, "-snv"),
                RefLen > AltLen ~ paste0(RefLen-1, "-del"),
                RefLen < AltLen ~ paste0(AltLen-1, "-ins"),
              )
            )
            df_mutation = df_mutation %>% dplyr::select(-RefLen, -AltLen)
            
            fun_genome <- function(x,y,z){
              syscom = paste0("samtools faidx -n 100 /tmp/MicroSEC/reference/Homo_sapiens_assembly38.fasta '",x,":",y,"-",z,"' - ")
              r <- system(syscom,ignore.stderr =T,intern = T)
              return(r)
            }
            
            df_mutation$mut_len = as.integer(str_split(df_mutation$Mut_type, "-",simplify=TRUE)[,1])
            df_mutation$mut_type = str_split(df_mutation$Mut_type, "-",simplify=TRUE)[,2]
            
            df_mutation = df_mutation %>% dplyr::mutate(
              Neighbor_start_1 = case_when(
                mut_type == "ins" ~ Pos - 19,
                mut_type == "del" ~ Pos - 19,
                mut_type == "snv" ~ Pos - 20
              ),
              Neighbor_end_1 = Pos - 1,
              Neighbor_start_2 = case_when(
                mut_type == "ins" ~ Pos + 1,
                mut_type == "del" ~ Pos + mut_len + 1,
                mut_type == "snv" ~ Pos + mut_len
              ),
              Neighbor_end_2 = case_when(
                mut_type == "ins" ~ Pos + 20,
                mut_type == "del" ~ Pos + mut_len + 20,
                mut_type == "snv" ~ Pos + mut_len + 19
              )
            )
            df_mutation$Pre_Neighbor = ""
            df_mutation$Post_Neighbor = ""
            for(i in seq_len(length(df_mutation$Chr))){
              df_mutation$Pre_Neighbor[i] = fun_genome(df_mutation$Chr[i],df_mutation$Neighbor_start_1[i],df_mutation$Neighbor_end_1[i])[2]
              df_mutation$Post_Neighbor[i] = fun_genome(df_mutation$Chr[i],df_mutation$Neighbor_start_2[i],df_mutation$Neighbor_end_2[i])[2]
            }
            df_mutation = df_mutation %>% dplyr::mutate(
              Neighborhood_sequence = paste0(Pre_Neighbor, Alt, Post_Neighbor)
            )
            df_mutation = df_mutation %>% dplyr::select(-Pre_Neighbor, -Post_Neighbor)
            df_mutation = df_mutation[order(df_mutation$Chr, df_mutation$Pos),]
            df_mutation = df_mutation %>% dplyr::select(-mut_len, -mut_type,-Neighbor_start_1, -Neighbor_end_1, -Neighbor_start_2, -Neighbor_end_2)
            write_tsv(progress = F,df_mutation,mutation_file)
          }
        }
        
        print("Trimming BAM/CRAM file...")
        if(!file.exists(output) & dim(tmp)[1]!=0){
          rm(tmp)
          if(!file.exists(paste0(bam_file_slim, ".bai"))){
            df_mutation_save = df_mutation
            download_region = data.frame(matrix(rep(NA,3),nrow=1))[numeric(0),]
            colnames(download_region) = c("chrom", "chromStart", "chromEnd")
            print(paste0(sample_name, ": ", dim(df_mutation_save)[1], " mutations"))
            for(i in chromosomes){
              df_mutation = df_mutation_save[df_mutation_save$Chr == i,]
              continuous = FALSE
              for(mut_no in seq_len(dim(df_mutation)[1])){
                if(mut_no == 1 & mut_no != dim(df_mutation)[1]){
                  if(df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
                     df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400){
                    download_region = rbind(download_region, c(df_mutation$Chr[mut_no],
                                                               max(1, df_mutation$Pos[mut_no] - 200) - 1,
                                                               df_mutation$Pos[mut_no] + 200))
                    colnames(download_region) = c("chrom", "chromStart", "chromEnd")
                    continuous = FALSE
                  } else{
                    continuous = TRUE
                    pos_last = max(1, df_mutation$Pos[mut_no] - 200)
                  }
                } else if (mut_no == 1 & mut_no == dim(df_mutation)[1]){
                  download_region = rbind(download_region, c(df_mutation$Chr[mut_no],
                                                             max(1, df_mutation$Pos[mut_no] - 200) - 1,
                                                             df_mutation$Pos[mut_no] + 200))
                  colnames(download_region) = c("chrom", "chromStart", "chromEnd")
                } else if (mut_no == dim(df_mutation)[1]){
                  if(continuous){
                    download_region = rbind(download_region, c(df_mutation$Chr[mut_no],
                                                               pos_last - 1,
                                                               df_mutation$Pos[mut_no] + 200))
                    colnames(download_region) = c("chrom", "chromStart", "chromEnd")
                  } else{
                    download_region = rbind(download_region, c(df_mutation$Chr[mut_no],
                                                               max(1, df_mutation$Pos[mut_no] - 200) - 1,
                                                               df_mutation$Pos[mut_no] + 200))
                    colnames(download_region) = c("chrom", "chromStart", "chromEnd")
                  }
                } else{
                  if(continuous){
                    if(df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
                       df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400){
                      download_region = rbind(download_region, c(df_mutation$Chr[mut_no],
                                                                 pos_last - 1,
                                                                 df_mutation$Pos[mut_no] + 200))
                      colnames(download_region) = c("chrom", "chromStart", "chromEnd")
                      continuous = FALSE
                    }
                  } else{
                    if(df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
                       df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400){
                      download_region = rbind(download_region, c(df_mutation$Chr[mut_no],
                                                                 max(1, df_mutation$Pos[mut_no] - 200) - 1,
                                                                 df_mutation$Pos[mut_no] + 200))
                      colnames(download_region) = c("chrom", "chromStart", "chromEnd")
                    } else{
                      continuous = TRUE
                      pos_last = max(1, df_mutation$Pos[mut_no] - 200)
                    }
                  }
                }
              }
            }
            write_tsv(progress = F,col_names = F,download_region, paste0(bam_file,".bed"))
            rm(download_region)
            system_out = 1
            syscom = paste0("samtools view -h --no-PG ",
                            bam_file, " -T ",
                            #url_cram, " -T ",
                            reference_genome, " ",
                            "-X ", file_crai, " -L ",
                            paste0(bam_file,".bed"), " > ",
                            bam_file_slim)
            system_out = (system(syscom))
            if(system_out == 0){
              system_out = 1
              syscom = paste0("samtools sort -o ", bam_file_tmp, " ", bam_file_slim)
              print("Sorting BAM file...")
              system_out = (system(syscom))
              if(system_out == 0){
                system_out = 1
                syscom = paste0("samtools view -bS ", bam_file_tmp, " > ", bam_file_slim)
                print("Compressing BAM file...")
                system_out = (system(syscom))
                if(system_out == 0){
                  system_out = 1
                  syscom = paste0("samtools index ", bam_file_slim)
                  print("Indexing BAM file...")
                  system_out = (system(syscom))
                  if(system_out == 0){
                    print(paste0("Slimmed BAM files were saved as ", bam_file_slim))
                    file.remove(bam_file)
                    file.remove(bam_file_tmp)
                  }
                }
              }
            }
            df_mutation = df_mutation_save
            rm(df_mutation_save)
          }
          
          print("FFPE artifact filtering by MicroSEC...")
          if(!file.exists(output) & file.exists(paste0(bam_file_slim, ".bai"))){
            bam_file = bam_file_slim
            
            fun_load_bam(bam_file)
            
            result = fun_read_check(short_homology_search_length = 4)
            msec = result[[1]]
            homology_search = result[[2]]
            mut_depth = result[[3]]
            rm(result)
            msec_tmp = fun_homology(msec,homology_search, min_homology_search = 40)
            msec_tmp = fun_summary(msec_tmp)
            msec_tmp = fun_analysis(msec_tmp, mut_depth,
                                         short_homology_search_length = 4,
                                         min_homology_search = 40,
                                         threshold_p = 10 ^ (-6),
                                         threshold_hairpin_ratio = 0.50,
                                         threshold_short_length = 0.75,
                                         threshold_distant_homology = 0.15,
                                         threshold_low_quality_rate = 0.1,
                                         homopolymer_length = 15)
            
            distant_homology_rate_40 = msec_tmp$distant_homology_rate
            not_long_repeat_40 = msec_tmp$not_long_repeat
            
            msec_tmp = fun_homology(msec,homology_search, min_homology_search = 35)
            msec_tmp = fun_summary(msec_tmp)
            msec_tmp = fun_analysis(msec_tmp, mut_depth,
                                         short_homology_search_length = 4,
                                         min_homology_search = 35,
                                         threshold_p = 10 ^ (-6),
                                         threshold_hairpin_ratio = 0.50,
                                         threshold_short_length = 0.75,
                                         threshold_distant_homology = 0.15,
                                         threshold_low_quality_rate = 0.1,
                                         homopolymer_length = 15)
            
            distant_homology_rate_35 = msec_tmp$distant_homology_rate
            not_long_repeat_35 = msec_tmp$not_long_repeat
            
            msec_tmp = fun_homology(msec,homology_search, min_homology_search = 30)
            msec_tmp = fun_summary(msec_tmp)
            msec_tmp = fun_analysis(msec_tmp, mut_depth,
                                         short_homology_search_length = 4,
                                         min_homology_search = 30,
                                         threshold_p = 10 ^ (-6),
                                         threshold_hairpin_ratio = 0.50,
                                         threshold_short_length = 0.75,
                                         threshold_distant_homology = 0.15,
                                         threshold_low_quality_rate = 0.1,
                                         homopolymer_length = 15)
            
            distant_homology_rate_30 = msec_tmp$distant_homology_rate
            not_long_repeat_30 = msec_tmp$not_long_repeat
            
            msec = fun_homology(msec,homology_search, min_homology_search = 25)
            msec = fun_summary(msec)
            msec = fun_analysis(msec, mut_depth,
                                     short_homology_search_length = 4,
                                     min_homology_search = 25,
                                     threshold_p = 10 ^ (-6),
                                     threshold_hairpin_ratio = 0.50,
                                     threshold_short_length = 0.75,
                                     threshold_distant_homology = 0.15,
                                     threshold_low_quality_rate = 0.1,
                                     homopolymer_length = 15)
            
            distant_homology_rate_25 = msec$distant_homology_rate
            not_long_repeat_25 = msec$not_long_repeat
            
            msec$distant_homology_rate_40 = distant_homology_rate_40
            msec$not_long_repeat_40 = not_long_repeat_40
            msec$distant_homology_rate_35 = distant_homology_rate_35
            msec$not_long_repeat_35 = not_long_repeat_35
            msec$distant_homology_rate_30 = distant_homology_rate_30
            msec$not_long_repeat_30 = not_long_repeat_30
            msec$distant_homology_rate_25 = distant_homology_rate_25
            msec$not_long_repeat_25 = not_long_repeat_25
            
            fun_save_gz(msec, output)
            rm(msec, msec_tmp, homology_search, mut_depth,
               distant_homology_rate_25, distant_homology_rate_30,
               distant_homology_rate_35, distant_homology_rate_40,
               not_long_repeat_25, not_long_repeat_30,
               not_long_repeat_35, not_long_repeat_40,
               df_bam)
          }
        }
      }
    }
    file.remove(sample_name)
  }
}
