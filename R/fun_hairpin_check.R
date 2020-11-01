#' Hairpin-structure sequence check function
#'
#' This function attempts to find hairpin structure sequences.
#'
#' @param Hairpin_seq_tmp The sequence to be checked.
#' @param Ref_seq Reference sequence around the mutation.
#' @param Hairpin_length The temporal length of hairpin sequences.
#' @param hair The length of sequences to be checked.
#' @return list(Hairpin_length, whether hairpin sequences exist or not)
#' @importFrom Biostrings matchPattern
#' @importFrom Biostrings DNAString
#' @examples
#' fun_hairpin_check(DNAString("GGGGGGGGAAAAAAAA"), DNAString("GGGGGGGGAAAAAAAATATATATATTTTTTTTCCCCCCCC"), 15, 16)
#' @export
fun_hairpin_check = function(Hairpin_seq_tmp, Ref_seq, Hairpin_length, hair){
  Reverse_seq = matchPattern(pattern = Hairpin_seq_tmp, subject = Ref_seq, max.mismatch=0, min.mismatch=0,  with.indels=FALSE, fixed=TRUE)
  if(length(Reverse_seq) > 0){
    Hairpin_length = max(Hairpin_length, hair)
    return(list(Hairpin_length, 1))
  }
  else{
    return(list(Hairpin_length, 0))
  }
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
