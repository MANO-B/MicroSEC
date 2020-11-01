#' Mutated position search function.
#'
#' This function attempts to find the mutated bases in each read.
#'
#' @param PRE The 5' side bases of the sequence for searching.
#' @param POST The 3' side bases of the sequence for searching.
#' @param Neighbor_seq Short reference sequence around the mutation.
#' @param neighbor_length The length from the mutation to the ends of the short reference sequence.
#' @param Alt_length The length of altered bases.
#' @return list(Pre_search_length, Post_search_length, Peri_seq_1, Peri_seq_2)
#' @examples
#' fun_setting(0, 20, "AGATCCT", 3, 1)
#' @export
fun_setting = function(PRE, POST, Neighbor_seq, neighbor_length, Alt_length){
  Peri_seq_1 = substring(Neighbor_seq, neighbor_length + 1 - PRE, neighbor_length + POST + Alt_length)
  Peri_seq_2 = substring(Neighbor_seq, neighbor_length + 1 - POST, neighbor_length + PRE + Alt_length)
  return(list(PRE, POST, Peri_seq_1, Peri_seq_2))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
