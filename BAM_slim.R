# to reduce the size of BAM files by deleting parts not related to mutations.
# It is recommended to run this program when the size of the BAM file is more than 10GB.
# dependency: samtools
# how to use
# Rscript BAM_slim.R [sample information tsv file] [progress bar Y/N]
# load necessary packages
library(MicroSEC)

# set arguments
args <- commandArgs(trailingOnly = T)
sample_list <- args[1]

sample_info <- read.csv(sample_list,
                        header = FALSE,
                        stringsAsFactors = FALSE,
                        sep = "\t")

for (sample in seq_len(dim(sample_info)[1])) {
  sample_name = sample_info[sample, 1]
  mutation_file = sample_info[sample, 2]
  bam_file = sample_info[sample, 3]
  bam_file_bai = paste(bam_file, ".bai", sep="")
  print(paste(sample_name, "/", sample, "/", dim(sample_info)[1]))
  
  if (!file.exists(bam_file_bai)) {
    print("sorting BAM file")
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
  
  # load mutation information
  fun_load_mutation(mutation_file, sample_name) # df_mutation
  df_mutation = df_mutation[order(df_mutation$Chr, df_mutation$Pos),]
  sep_new = TRUE
  continuous = FALSE
  chr_last = ""
  pos_last = 0
  bam_file_tmp1 = paste(bam_file, ".tmp1", sep="")
  bam_file_tmp2 = paste(bam_file, ".tmp2", sep="")
  bam_file_slim = paste(bam_file, ".SLIM", sep="")
  for (mut_no in seq_len(dim(df_mutation)[1])) {
    print(paste(mut_no, "/", dim(df_mutation)[1]))
    if (mut_no == 1 & mut_no != dim(df_mutation)[1]) {
      if (df_mutation$Chr[mut_no + 1] != df_mutation$Chr[mut_no] |
          df_mutation$Pos[mut_no + 1] > df_mutation$Pos[mut_no] + 400) {
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
                       df_mutation$Chr[mut_no],
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
                       df_mutation$Chr[mut_no],
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
                         df_mutation$Chr[mut_no],
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
                         df_mutation$Chr[mut_no],
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
                         df_mutation$Chr[mut_no],
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
  syscom = paste("samtools sort -@ 4 -o ",
                 bam_file_tmp2,
                 " ",
                 bam_file_slim,
                 sep="")
  system(syscom)
  syscom = paste("samtools view -bS ",
                 bam_file_tmp2,
                 " > ",
                 bam_file_slim,
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
