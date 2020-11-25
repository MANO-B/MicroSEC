#' Mutated position search function.
#'
#' This function attempts to find the mutated bases in each read.
#'
#' @param pre The 5' side bases of the sequence for searching.
#' @param post The 3' side bases of the sequence for searching.
#' @param neighbor_seq Short reference sequence around the mutation.
#' @param neighbor_length The length from the mutation to the ends of the short
#'   reference sequence.
#' @param alt_length The length of altered bases.
#' @return list(pre_search_length, post_search_length, peri_seq_1, peri_seq_2)
fun_setting <- function(pre, post, neighbor_seq, neighbor_length, alt_length) {
  peri_seq_1 <- substring(neighbor_seq,
                         neighbor_length + 1 - pre,
                         neighbor_length + post + alt_length)
  peri_seq_2 <- substring(neighbor_seq,
                         neighbor_length + 1 - post,
                         neighbor_length + pre + alt_length)
  return(list(pre, post, peri_seq_1, peri_seq_2))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
