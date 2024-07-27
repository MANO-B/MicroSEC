#' Homology check function.
#'
#' This function attempts to search the homologous regions.
#'
#' @param msec Mutation filtering information.
#' @param df_distant Sequences to be checked.
#' @param min_homology_search Minimum length to define "homologous".
#' @param ref_genome Reference genome for the data.
#' @param chr_no Reference genome chromosome number (human=24, mouse=22).
#' @param progress_bar "Y": You can see the progress visually.
#' @return msec
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom stringr str_sub
#' @importFrom stringr str_detect
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings PDict
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings countPDict
#' @importFrom BiocGenerics as.data.frame
#' @importFrom BiocGenerics paste
#' @importFrom GenomeInfoDb seqnames
#' @examples
#' \dontrun{
#' data(msec_read_checked)
#' data(homology_searched)
#' fun_homology(msec = msec_read_checked,
#'    df_distant = homology_searched,
#'    min_homology_search = 40,
#'    ref_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#'    chr_no = 24,
#'    progress_bar = "N"
#' )
#' }
#' @export
fun_homology <- function(msec,
                         df_distant,
                         min_homology_search,
                         ref_genome,
                         chr_no,
                         progress_bar) {
  # initialize
  Seq <- NULL
  Direction <- NULL
  fixed_seq_pre <- NULL
  fixed_seq_post <- NULL
  
  if (!is.null(df_distant)) {
    msec$distant_homology <- 0
    df_distant <- df_distant %>% dplyr::mutate(
      Seq = as.character(Seq))
    df_distant <- df_distant %>% dplyr::mutate(
      length = nchar(Seq))
    df_distant$number <- floor(seq(1, (dim(df_distant)[1] / 2 + 0.5), 0.5))
    max_no <- max(df_distant$number)
    df_distant <- df_distant[df_distant$length >= min_homology_search, ]
    df_distant <- df_distant %>% dplyr::mutate(
      fixed_seq_pre = str_sub(
        Seq, length - min_homology_search + 1, length),
      fixed_seq_post = str_sub(
        Seq, 1, min_homology_search)
    )
    df_distant <- df_distant %>% dplyr::mutate(
      fixed_seq = ifelse(Direction == "pre", fixed_seq_pre, fixed_seq_post)
    )
    df_distant <- df_distant[!str_detect(df_distant$fixed_seq, pattern = "N"), ]
    search_seq_f <- DNAStringSet(df_distant$fixed_seq)
    search_seq_r <- reverseComplement(search_seq_f)
    if (length(search_seq_f) > 0) {
      search_seq_f <- PDict(search_seq_f)
      search_seq_r <- PDict(search_seq_r)
      distant_homology <- rep(0, length(search_seq_f))
      if (progress_bar == "Y") {
        pb <- utils::txtProgressBar(min = 0,
                            max = chr_no,
                            width = 20,
                            style = 3)
      }
      for (seqname in 1:chr_no) {
        if (progress_bar == "Y") {
          utils::setTxtProgressBar(pb, seqname)
          cat(paste0("Chromosome screening: ", seqname, "/", chr_no, "    "))
        }
        target <- ref_genome[[seqnames(ref_genome)[[seqname]]]]
        distant_homology <- distant_homology + countPDict(search_seq_f, target)
        distant_homology <- distant_homology + countPDict(search_seq_r, target)
        gc()
        gc()
      }
      df_distant$distant_homology <- distant_homology
      if (progress_bar == "Y") {
        pb <- utils::txtProgressBar(min = 0,
                            max = max(1, max_no),
                            width = 20,
                            style = 3)
        pb_t <- ceiling(max_no / 100)
      }
      tmp_distant = df_distant[1,]
      tmp_distant$sample_name = ""
      tmp_distant$Chr = ""
      tmp_distant$Pos = 0
      tmp_distant$Ref = ""
      tmp_distant$Alt = ""
      k = 0
      if (dim(df_distant[df_distant$number == 1, ])[1] == 1) {
        tmp_distant <- df_distant[df_distant$number == 1, ][1, ]
        if (tmp_distant$distant_homology > 0) {
          k = 1
        } else {
          k = 0
        }
      }
      if (dim(df_distant[df_distant$number == 1, ])[1] == 2) {
        tmp_distant <- df_distant[df_distant$number == 1, ][1, ]
        tmp_distant_2 <- df_distant[df_distant$number == 1, ][2, ]
        if (tmp_distant$distant_homology > 0 |
            tmp_distant_2$distant_homology > 0) {
          k = 1
        } else {
          k = 0
        }
      }
      Sample_previous = tmp_distant$sample_name
      Chr_previous = tmp_distant$Chr
      Pos_previous = tmp_distant$Pos
      Ref_previous = tmp_distant$Ref
      Alt_previous = tmp_distant$Alt
      if (max_no > 1){
        for (i in 2:max_no) {
          if (progress_bar == "Y") {
            utils::setTxtProgressBar(pb, i)
            if ((i - 1) %% pb_t == 0) {
              cat(paste0("Homology count: ", i, "/", max_no, "      "))
            }
          }
          if (dim(df_distant[df_distant$number == i, ])[1] == 1) {
            tmp_distant <- df_distant[df_distant$number == i, ][1, ]
            if (tmp_distant$sample_name != Sample_previous |
                tmp_distant$Chr != Chr_previous |
                tmp_distant$Pos != Pos_previous |
                tmp_distant$Ref != Ref_previous |
                tmp_distant$Alt != Alt_previous) {
              if (k > 0) {
                index_no = which(msec$Sample == Sample_previous &
                                        msec$Chr == Chr_previous &
                                        msec$Pos == Pos_previous &
                                        msec$Ref == Ref_previous &
                                        msec$Alt == Alt_previous)
                msec[index_no, ]$distant_homology <-
                  msec[index_no, ]$distant_homology + k
              }
              Sample_previous = tmp_distant$sample_name
              Chr_previous = tmp_distant$Chr
              Pos_previous = tmp_distant$Pos
              Ref_previous = tmp_distant$Ref
              Alt_previous = tmp_distant$Alt
              k = 0
            }
            if (tmp_distant$distant_homology > 0) {
              k = k + 1
            }
          }
          if (dim(df_distant[df_distant$number == i, ])[1] == 2) {
            tmp_distant <- df_distant[df_distant$number == i, ][1, ]
            tmp_distant_2 <- df_distant[df_distant$number == i, ][2, ]
            if (tmp_distant$sample_name != Sample_previous |
                tmp_distant$Chr != Chr_previous |
                tmp_distant$Pos != Pos_previous |
                tmp_distant$Ref != Ref_previous |
                tmp_distant$Alt != Alt_previous) {
              if (k > 0) {
                index_no = which(msec$Sample == Sample_previous &
                                   msec$Chr == Chr_previous &
                                   msec$Pos == Pos_previous &
                                   msec$Ref == Ref_previous &
                                   msec$Alt == Alt_previous)
                msec[index_no, ]$distant_homology <-
                  msec[index_no, ]$distant_homology + k
              }
              Sample_previous = tmp_distant$sample_name
              Chr_previous = tmp_distant$Chr
              Pos_previous = tmp_distant$Pos
              Ref_previous = tmp_distant$Ref
              Alt_previous = tmp_distant$Alt
              k = 0
            }
            if (tmp_distant$distant_homology > 0 |
                tmp_distant_2$distant_homology > 0) {
              k = k + 1
            }
          }
        }  
      }
      if (k > 0) {
        index_no = which(msec$Sample == Sample_previous &
                           msec$Chr == Chr_previous &
                           msec$Pos == Pos_previous &
                           msec$Ref == Ref_previous &
                           msec$Alt == Alt_previous)
        msec[index_no, ]$distant_homology <-
          msec[index_no, ]$distant_homology + k
      }
    }
  }  
  return(msec)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
