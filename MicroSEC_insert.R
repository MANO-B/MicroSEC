# Insert length checker

# load necessary packages
library(MicroSEC)

# set arguments
args = commandArgs(trailingOnly = T)

if(args[3] == "N" | args[3] == "Y"){
  wd = args[1]
  SAMPLE_LIST = args[2]

  setwd(wd)
  
  # load sample information tsv file
  SAMPLE_INFO = read.csv(SAMPLE_LIST,
                         header=FALSE,
                         stringsAsFactors=FALSE,
                         sep="\t")
  
  # initialize
  INSERT = NULL
    
  for(SAMPLE in 1:dim(SAMPLE_INFO)[1]){
    SAMPLE_NAME = SAMPLE_INFO[SAMPLE,1]
    BAM_FILE = SAMPLE_INFO[SAMPLE,3]

    # load mutation information
    insert_length = fun_insert_length(BAM_FILE)
    tmp = data.frame(SAMPLE_NAME, insert_length)
    colnames(tmp) = c(Sample, Insert_length)
    INSERT = rbind(INSERT, tmp)
  }
  # save the results
  write.table(INSERT, file=paste(wd, "/Insert_size_", SAMPLE_INFO[1,1], ".tsv", sep=""), sep = "\t", na="", row.names=FALSE, col.names=FALSE, quote=FALSE)
}else {
  SAMPLE_NAME = args[2]
  BAM_FILE = args[4]

  insert_length = fun_insert_length(BAM_FILE)
  tmp = data.frame(SAMPLE_NAME, insert_length)
  colnames(tmp) = c(Sample, Insert_length)

  # save the results
  write.table(tmp, file=paste("Insert_size_", SAMPLE_NAME, ".tsv", sep=""), sep = "\t", na="", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

