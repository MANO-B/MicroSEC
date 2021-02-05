#' Supporting length calculation function.
#'
#' This function attempts to calculate supporting lengths of a read.
#'
#' @param df_cigar The CIGAR data of the read.
#' @param df_pos The start position of the read.
#' @param mut_position The mutation position in the genome.
#' @param alt_length The length of altered bases.
#' @param indel_status The mutation is indel or not.
#' @importFrom stringr str_split
#' @importFrom BiocGenerics as.data.frame
#' @return list(pre_support_length, post_support_length, short_support_length)
fun_normal_support <- function(df_cigar,
                               df_pos,
                               mut_position,
                               alt_length,
                               indel_status) {
  cigar_num <- as.integer(str_split(df_cigar, "[:upper:]")[[1]])
  cigar_type <- str_split(df_cigar, "[:digit:]+")[[1]][-1]
  tmp_pos <- df_pos - 1
  alt <- alt_length * (1 - indel_status)
  for (k in seq_len(length(cigar_type))) {
    if (cigar_type[k] == "D") {
      tmp_pos <- tmp_pos + cigar_num[k]
    }
    if (cigar_type[k] == "M") {
      tmp_pos <- tmp_pos + cigar_num[k]
    }
  }
  if (tmp_pos < (mut_position + alt)) {
    return(list(0, 0, -1))
  } else {
    return(list(min(mut_position - df_pos, 199),
                max(0, min(tmp_pos - mut_position, 199)),
                max(0, min(mut_position - df_pos, tmp_pos - mut_position, 99))))
  }
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
