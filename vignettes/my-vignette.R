## ----setting------------------------------------------------------------------
wd <- "~"
knitr::opts_chunk$set(collapse = TRUE,
                      fig.width = 12,
                      fig.height = 8,
                      echo = TRUE,
                      warning = FALSE,
                      message = TRUE,
                      comment = "#>")
options(rmarkdown.html_vignette.check_title = FALSE,
        show.error.messages = FALSE,
        warn = -1)

progress_bar <- "N"

## ----packages-----------------------------------------------------------------
library(MicroSEC)

## ----analysis-----------------------------------------------------------------
# initialize
msec <- NULL
homology_search <- NULL
mut_depth <- NULL

# test data
sample_name <- "sample"
read_length <- 150
adapter_1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
organism <- "hg38"

# load mutation information
df_mutation <- fun_load_mutation(
   system.file("extdata", "mutation_list.tsv", package = "MicroSEC"),
   "sample",
   BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
   24)
df_bam <- fun_load_bam(
   system.file("extdata", "sample.bam", package = "MicroSEC"))

# another example data
# data(exampleMutation)
# data(exampleBam)
# df_mutation <- exampleMutation
# df_bam <- exampleBam

# load genomic sequence
ref_genome <- fun_load_genome(organism)
chr_no <- fun_load_chr_no(organism)

# analysis
result <- fun_read_check(df_mutation = df_mutation,
                         df_bam = df_bam,
                         ref_genome = ref_genome,
                         sample_name = sample_name,
                         read_length = read_length,
                         adapter_1 = adapter_1,
                         adapter_2 = adapter_2,
                         short_homology_search_length = 4,
                         min_homology_search = 40,
                         progress_bar = progress_bar)
msec_read_checked <- result[[1]]
homology_searched <- result[[2]]
mut_depth_checked <- result[[3]]

# search homologous sequences
msec_homology_searched = fun_homology(msec = msec_read_checked,
                    df_distant = homology_searched,
                    min_homology_search = 40,
                    ref_genome = ref_genome,
                    chr_no = chr_no,
                    progress_bar = progress_bar)

# statistical analysis
msec_summarized <- fun_summary(msec_homology_searched)
msec_analyzed <- fun_analysis(msec = msec_summarized,
                    mut_depth = mut_depth_checked,
                    short_homology_search_length = 4,
                    min_homology_search = 40,
                    threshold_p = 10 ^ (-6),
                    threshold_hairpin_ratio = 0.50,
                    threshold_short_length = 0.75,
                    threshold_distant_homology = 0.15,
                    threshold_soft_clip_ratio = 0.50,
                    threshold_low_quality_rate = 0.1,
                    homopolymer_length = 15)

# save the results as a tsv.gz file.
#fun_save(msec_analyzed, "~/MicroSEC_test.tsv.gz")

## ----result-------------------------------------------------------------------
msec_analyzed

## ----sessioninfo--------------------------------------------------------------
sessionInfo()

