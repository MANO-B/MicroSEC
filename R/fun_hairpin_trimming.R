#' Adapter sequence trimming function
#'
#' This function attempts to find and cut the adapter sequences in the ends.
#'
#' @param hairpin_seq The sequence to be trimmed.
#' @param mut_read_strand The strand of the sequence, "+" or "-".
#' @param adapter_1 The Read 1 adapter sequence of the library.
#' @param adapter_2 The Read 2 adapter sequence of the library.
#' @return Adapter-trimmed hairpin_seq
#' @importFrom Biostrings trimLRPatterns
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom BiocGenerics as.data.frame
#' @importFrom stringr str_split
fun_hairpin_trimming <- function(hairpin_seq,
                                mut_read_strand,
                                adapter_1,
                                adapter_2) {
  if (mut_read_strand == "+") {
    hairpin_seq <- trimLRPatterns(Rpattern = adapter_1,
                                 subject = hairpin_seq,
                                 max.Rmismatch = 0.1)
    hairpin_seq <- trimLRPatterns(Rpattern = adapter_2,
                                 subject = hairpin_seq,
                                 max.Rmismatch = 0.1)
    hairpin_seq <- str_split(hairpin_seq,
                            str_sub(adapter_1, 1, 15))[[1]]
    if (length(hairpin_seq) > 1) {
      hairpin_seq <- hairpin_seq[[1]]
    }
    hairpin_seq <- str_split(hairpin_seq,
                            str_sub(adapter_2, 1, 15))[[1]]
    if (length(hairpin_seq) > 1) {
      hairpin_seq <- hairpin_seq[[1]]
    }
    hairpin_seq <- reverseComplement(DNAString(hairpin_seq))
  }
  else{
    hairpin_seq <- reverseComplement(hairpin_seq)
    hairpin_seq <- trimLRPatterns(Rpattern = adapter_1,
                                 subject = hairpin_seq,
                                 max.Rmismatch = 0.1)
    hairpin_seq <- trimLRPatterns(Rpattern = adapter_2,
                                 subject = hairpin_seq,
                                 max.Rmismatch = 0.1)
    hairpin_seq <- str_split(hairpin_seq,
                            str_sub(adapter_1, 1, 15))[[1]]
    if (length(hairpin_seq) > 1) {
      hairpin_seq <- hairpin_seq[[1]]
    }
    hairpin_seq <- str_split(hairpin_seq,
                            str_sub(adapter_2, 1, 15))[[1]]
    if (length(hairpin_seq) > 1) {
      hairpin_seq <- hairpin_seq[[1]]
    }
    hairpin_seq <- DNAString(hairpin_seq)
  }
  return(hairpin_seq)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
