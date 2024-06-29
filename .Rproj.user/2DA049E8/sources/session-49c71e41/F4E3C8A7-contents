#' Supporting length calculation function.
#'
#' This function attempts to calculate supporting lengths of a read.
#'
#' @param df_cigar The CIGAR data of the read.
#' @param df_seq The sequence to be checked.
#' @param mut_read_strand The strand of the sequence, "+" or "-".
#' @param adapter_1 The Read 1 adapter sequence of the library.
#' @param adapter_2 The Read 2 adapter sequence of the library.
#' @param mut_position The mutation position in the read.
#' @param alt_length The length of altered bases.
#' @param indel_status The mutation is indel or not.
#' @importFrom stringr str_split
#' @importFrom BiocGenerics as.data.frame
#' @return list(pre_support_length, post_support_length, soft_clipped_read)
fun_support <- function(df_cigar,
                       df_seq,
                       mut_read_strand,
                       adapter_1,
                       adapter_2,
                       mut_position,
                       alt_length,
                       indel_status) {
  pre_support_length_tmp <- 0
  post_support_length_tmp <- 0
  cigar_num <- as.integer(str_split(df_cigar, "[:upper:]")[[1]])
  cigar_type <- str_split(df_cigar, "[:digit:]+")[[1]][-1]
  tmp_pos <- 1
  soft_clipped_read_tmp <- 0
  for (k in seq_len(length(cigar_type))) {
    if (cigar_type[k] == "H") {
      tmp_pos <- tmp_pos + cigar_num[k]
    }
    if (cigar_type[k] == "I") {
      tmp_pos <- tmp_pos + cigar_num[k]
    }
    if (cigar_type[k] == "S") {
      clipped_seq <- fun_hairpin_trimming(
        df_seq[tmp_pos:(tmp_pos + cigar_num[k] - 1)],
        mut_read_strand,
        adapter_1, adapter_2)
      if (length(clipped_seq) > 0) {
        soft_clipped_read_tmp <- 1
      }
      tmp_pos <- tmp_pos + cigar_num[k]
    }
    if (cigar_type[k] == "M") {
      if ((tmp_pos + cigar_num[k] - 1) < mut_position) {
        pre_support_length_tmp <-
          pre_support_length_tmp + cigar_num[k]
      } else if (tmp_pos > (mut_position + alt_length - 1)) {
        post_support_length_tmp <-
          post_support_length_tmp + cigar_num[k]
      } else {
        if (indel_status == 1) {
          pre_support_length_tmp <-
            pre_support_length_tmp + (mut_position - tmp_pos)
          post_support_length_tmp <-
            post_support_length_tmp +
            (cigar_num[k] + tmp_pos - mut_position)
        } else {
          pre_support_length_tmp <-
            pre_support_length_tmp + (mut_position - tmp_pos)
          post_support_length_tmp <-
            post_support_length_tmp +
            (cigar_num[k] + tmp_pos - mut_position - alt_length)
        }
      }
      tmp_pos <- tmp_pos + cigar_num[k]
    }
  }
  return(list(pre_support_length_tmp,
              post_support_length_tmp,
              soft_clipped_read_tmp))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
