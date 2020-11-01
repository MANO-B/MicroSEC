## ---- eval=FALSE--------------------------------------------------------------
#  wd = "/mnt/HDD8TB/MicroSEC"
#  knitr::opts_chunk$set(collapse = TRUE, fig.width=12, fig.height=8, echo=TRUE, warning=FALSE, message=TRUE, comment = "#>")
#  knitr::opts_knit$set(root.dir = wd)
#  options(rmarkdown.html_vignette.check_title = FALSE)

## ----packages, eval=FALSE-----------------------------------------------------
#  library(tidyr)
#  library(openxlsx)
#  library(data.table)
#  library(R.utils)
#  library(stringr)
#  library(magrittr)
#  library(dplyr)
#  library(Biostrings)
#  library(BSgenome.Hsapiens.UCSC.hg38)
#  library(BSgenome.Mmusculus.UCSC.mm10)
#  library(GenomicAlignments)
#  library(Rsamtools)
#  library(tcltk)
#  library(gtools)
#  library(MicroSEC)

## ---- eval=FALSE--------------------------------------------------------------
#  # load sample information tsv file
#  SAMPLE_INFO = read.csv(paste(wd, "/source/Sample_list.txt", sep=""), header=FALSE, stringsAsFactors=FALSE, sep="\t")
#  
#  for(SAMPLE in 1:dim(SAMPLE_INFO)[1]){
#    SAMPLE_NAME = SAMPLE_INFO[SAMPLE,1]
#    MUTATION_FILE = SAMPLE_INFO[SAMPLE,2]
#    BAM_FILE = SAMPLE_INFO[SAMPLE,3]
#    MUTATION_SUPPORTING_READ_LIST = SAMPLE_INFO[SAMPLE,4]
#    READ_length = as.integer(SAMPLE_INFO[SAMPLE,5])
#    ADAPTOR_SEQ = SAMPLE_INFO[SAMPLE,6]
#    GENOME = SAMPLE_INFO[SAMPLE,7]
#  
#    # initialize
#    MSEC = NULL
#    Homology_search = NULL
#  
#    # load mutation information
#    df_mutation = fun_load_mutation(MUTATION_FILE)
#    df_BAM = fun_load_BAM(BAM_FILE)
#    df_mut_call = fun_load_ID(MUTATION_SUPPORTING_READ_LIST)
#  
#    # load genomic sequence
#    genome = fun_load_genome(GENOME)
#    Chr_No = fun_load_chr_no(GENOME)
#  
#    # analysis
#    result = fun_read_check(df_mutation = df_mutation,
#                            df_BAM =  df_BAM,
#                            df_mut_call = df_mut_call,
#                            genome = genome,
#                            Chr_No = Chr_No,
#                            SAMPLE_NAME = SAMPLE_NAME,
#                            READ_length = READ_length,
#                            ADAPTOR_SEQ = ADAPTOR_SEQ,
#                            PROGRESS_BAR = "Y")
#    MSEC = rbind(MSEC, result[[1]])
#    Homology_search = rbind(Homology_search, result[[2]])
#  
#    # search homologous sequences
#    MSEC = fun_homology(MSEC,
#                          Homology_search,
#                          PROGRESS_BAR = "Y")
#  
#  # statistical analysis
#    MSEC = fun_summary(MSEC)
#  
#  # save the results
#    fun_save(MSEC, wd)
#  }

## -----------------------------------------------------------------------------
sessionInfo()

