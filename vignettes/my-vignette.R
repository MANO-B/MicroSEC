## -----------------------------------------------------------------------------
wd = "."
knitr::opts_chunk$set(collapse = TRUE,
                      fig.width = 12,
                      fig.height = 8,
                      echo = TRUE,
                      warning = FALSE,
                      message = TRUE,
                      comment = "#>")
options(rmarkdown.html_vignette.check_title = FALSE)
progress_bar = "N"

## ----packages-----------------------------------------------------------------
library(MicroSEC)

## -----------------------------------------------------------------------------
# initialize
msec = NULL
homology_search = NULL
mut_depth = NULL

# test data
sample_name = "H15-11943-1-T_TDv3"
read_length = 151
adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
organism = "hg38"

# load mutation information
data(exampleMutation)
data(exampleBAM)
data(exampleMutCall)

# load genomic sequence
ref_genome = fun_load_genome(organism)
chr_no = fun_load_chr_no(organism)

# analysis
result = fun_read_check(df_mutation = exampleMutation,
                        df_bam = exampleBAM,
                        df_mut_call = exampleMutCall,
                        ref_genome = ref_genome,
                        sample_name = sample_name,
                        read_length = read_length,
                        adapter_1 = adapter_1,
                        adapter_2 = adapter_2,
                        short_homology_search_length = 4,
                        progress_bar = progress_bar)
msec = rbind(msec, result[[1]])
homology_search = rbind(homology_search, result[[2]])
mut_depth = rbind(mut_depth, result[[3]])

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
                    threshold_soft_clip_ratio = 0.90,
                    threshold_short_length = 0.8,
                    threshold_distant_homology = 0.2,
                    threshold_low_quality_rate = 0.1,
                    homopolymer_length = 15)

# save the results in the working/output directory.
# fun_save(msec, sample_name, wd)

## -----------------------------------------------------------------------------
msec

## -----------------------------------------------------------------------------
sessionInfo()

