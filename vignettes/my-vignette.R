## -----------------------------------------------------------------------------
wd = "/mnt/HDD8TB/MicroSEC/archive"
knitr::opts_chunk$set(collapse = TRUE, fig.width=12, fig.height=8, echo=TRUE, 
                      warning=FALSE, message=TRUE, comment = "#>")
knitr::opts_knit$set(root.dir = wd)
options(rmarkdown.html_vignette.check_title = FALSE)
PROGRESS_BAR = "N"

## ----packages-----------------------------------------------------------------
library(MicroSEC)

## -----------------------------------------------------------------------------
# initialize
MSEC = NULL
Homology_search = NULL
Mut_depth = NULL

# test data
SAMPLE_NAME = "H15-11943-1-T_TDv3"
READ_length = 151
ADAPTER_SEQ_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_SEQ_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
GENOME = "hg38"

# load mutation information
data(exampleMutation)
data(exampleBAM)
data(exampleMutCall)

# load genomic sequence
genome = fun_load_genome(GENOME)
Chr_No = fun_load_chr_no(GENOME)

# analysis
result = fun_read_check(df_mutation = exampleMutation,
                        df_BAM = exampleBAM,
                        df_mut_call = exampleMutCall,
                        genome = genome,
                        Chr_No = Chr_No,
                        SAMPLE_NAME = SAMPLE_NAME,
                        READ_length = READ_length,
                        ADAPTER_SEQ_1 = ADAPTER_SEQ_1,
                        ADAPTER_SEQ_2 = ADAPTER_SEQ_2,
                        Short_Homology_search_length = 4,
                        PROGRESS_BAR = PROGRESS_BAR)
MSEC = rbind(MSEC, result[[1]])
Homology_search = rbind(Homology_search, result[[2]])
Mut_depth = rbind(Mut_depth, result[[3]])

# search homologous sequences
MSEC = fun_homology(MSEC,
                    Homology_search,
                    Minimum_Homology_search_length = 40,
                    PROGRESS_BAR = PROGRESS_BAR)

# statistical analysis
MSEC = fun_summary(MSEC)
MSEC = fun_analysis(MSEC,
                    Mut_depth,
                    Short_Homology_search_length = 4,
                    Minimum_Homology_search_length = 40,
                    threshold_p = 10^(-6),
                    threshold_hairpin_ratio = 0.50,
                    threshold_soft_clip_ratio = 0.90,
                    threshold_short_length = 0.8,
                    threshold_distant_homology = 0.2,
                    threshold_low_quality_rate = 0.1,
                    Homopolymer_length = 15)

# save the results
fun_save(MSEC, SAMPLE_NAME, wd)

## -----------------------------------------------------------------------------
MSEC

## -----------------------------------------------------------------------------
sessionInfo()

