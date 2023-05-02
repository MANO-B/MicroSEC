#' Hairpin-structure sequence check function
#'
#' This function attempts to find hairpin structure sequences.
#'
#' @param hairpin_seq_tmp The sequence to be checked.
#' @param ref_seq Reference sequence around the mutation.
#' @param hairpin_length The temporal length of hairpin sequences.
#' @param hair The length of sequences to be checked.
#' @return list(hairpin_length, whether hairpin sequences exist or not)
#' @importFrom Biostrings matchPattern
#' @importFrom Biostrings DNAString
#' @importFrom BiocGenerics as.data.frame
fun_hairpin_check <- function(hairpin_seq_tmp, ref_seq, hairpin_length, hair) {
  reverse_seq <- matchPattern(pattern = hairpin_seq_tmp,
                             subject = ref_seq,
                             max.mismatch = 0,
                             min.mismatch = 0,
                             with.indels = FALSE,
                             fixed = TRUE)
  if (length(reverse_seq) > 0) {
    hairpin_length <- max(hairpin_length, hair)
    return(list(hairpin_length, 1))
  }
  else{
    return(list(hairpin_length, 0))
  }
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
