# Insert length checker

# load necessary packages
library(MicroSEC)

# set arguments
args = commandArgs(trailingOnly = T)

if(args[3] == "N" | args[3] == "Y"){
  wd = args[1]
  sample_list = args[2]

  setwd(wd)
  
  # load sample information tsv file
  sample_info = read.csv(sample_list,
                         header=FALSE,
                         stringsAsFactors=FALSE,
                         sep="\t")
  
  # initialize
  insert = NULL
    
  for(sample in 1:dim(sample_info)[1]){
    sample_name = sample_info[sample,1]
    bam_file = sample_info[sample,3]

    # load mutation information
    insert_length = fun_insert_length(bam_file)
    tmp = data.frame(sample_name, insert_length)
    colnames(tmp) = c("Sample", "Insert_length")
    insert = rbind(insert, tmp)
  }
  # save the results
  write.table(insert, file=paste(wd, "/Insert_size_", sample_info[1,1], ".tsv", sep=""), sep = "\t", na="", row.names=FALSE, col.names=FALSE, quote=FALSE)
}else {
  sample_name = args[2]
  bam_file = args[4]

  insert_length = fun_insert_length(bam_file)
  tmp = data.frame(sample_name, insert_length)
  colnames(tmp) = c("Sample", "Insert_length")

  # save the results
  write.table(tmp, file=paste("Insert_size_", sample_name, ".tsv", sep=""), sep = "\t", na="", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

