#' Adapter sequence trimming function
#'
#' This function attempts to find and cut the adapter sequences in the ends.
#'
#' @param Hairpin_seq The sequence to be trimmed.
#' @param mut_read_strand The strand of the sequence, "+" or "-".
#' @param ADAPTOR_SEQ The adapter sequence of the library.
#' @return Adapter-trimmed Hairpin_seq
#' @importFrom Biostrings trimLRPatterns
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom BiocGenerics as.data.frame
#' @importFrom stringr str_split
#' @examples
#' fun_hairpin_trimming(DNAString("GGAAAAAAATCTCTCAACA"), "+", "AGATCC")
#' @export
fun_hairpin_trimming = function(Hairpin_seq, mut_read_strand, ADAPTOR_SEQ){
  if(mut_read_strand == "+"){
    Hairpin_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ, subject = Hairpin_seq)
    Hairpin_seq = str_split(Hairpin_seq, ADAPTOR_SEQ)[[1]]
    if(length(Hairpin_seq)> 1){
      Hairpin_seq = Hairpin_seq[[1]]
    }
    Hairpin_seq = reverseComplement(DNAString(Hairpin_seq))
  }
  else{
    Hairpin_seq = reverseComplement(Hairpin_seq)
    Hairpin_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ, subject = Hairpin_seq)
    Hairpin_seq = str_split(Hairpin_seq, ADAPTOR_SEQ)[[1]]
    if(length(Hairpin_seq)> 1){
      Hairpin_seq = Hairpin_seq[[1]]
    }
    Hairpin_seq = DNAString(Hairpin_seq)
  }
  return(Hairpin_seq)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
