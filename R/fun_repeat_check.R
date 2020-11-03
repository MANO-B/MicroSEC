#' Repeat check function.
#'
#' This function attempts to check the repetitive sequence around the mutation.
#'
#' @param Rep_A The shorter sequence of [Ref] and [Alt].
#' @param Rep_B The longer sequence of [Ref] and [Alt].
#' @param Ref_seq Reference sequence around the mutation.
#' @param Width Search length for Ref_seq.
#' @param Del Insertion: 0, Deletion: 1
#' @return list(Pre_rep_status, Post_rep_status, Homopolymer_status)
#' @importFrom Biostrings DNAString
#' @importFrom BiocGenerics as.data.frame
#' @examples
#' fun_repeat_check(DNAString("A"), DNAString("AATC"), DNAString("GGAAAAAAATCTCTCAACA"), 5, 0)
#' @export
fun_repeat_check = function(Rep_A, Rep_B, Ref_seq, Width, Del){
  Pre_rep_status = 0
  Post_rep_status = 0
  Homopolymer_status = 0
  Rep_B = Rep_B[2:length(Rep_B)]
  homo_tmp_1 = 0
  homo_tmp_2 = 0
  for(q in 1:length(Rep_B)){
    Rep_seq = Rep_B[q:length(Rep_B)]
    check_rep = TRUE
    Post_rep_status_tmp = 0
    for(r in (1 + Del * length(Rep_B)):(Width - 1)){
      if(check_rep &  Ref_seq[Width + r + 1] == Rep_seq[((r - 1) %% length(Rep_seq) + 1)]){
        Post_rep_status_tmp = Post_rep_status_tmp + 1
      }
      else{
        check_rep = FALSE
      }
    }
    if(Post_rep_status_tmp >= length(Rep_seq)){
      Post_rep_status = max(Post_rep_status, Post_rep_status_tmp)
      if(gsub(as.character(Rep_seq[1]), "", as.character(Rep_seq)) == ""){
        homo_tmp_1 = Post_rep_status_tmp + length(Rep_seq)
      }
    }
  }
  for(q in 1:length(Rep_B)){
    Rep_seq = Rep_B[1:q]
    check_rep = TRUE
    Pre_rep_status_tmp = 0
    for(r in 1:(Width + 1)){
      if(check_rep &  Ref_seq[Width + 2 - r] == Rep_seq[((length(Rep_seq) - r) %% length(Rep_seq) + 1)]){
        Pre_rep_status_tmp = Pre_rep_status_tmp + 1
      }
      else{
        check_rep = FALSE
      }
    }
    if(Pre_rep_status_tmp >= length(Rep_seq)){
      Pre_rep_status = max(Pre_rep_status, Pre_rep_status_tmp)
      if(gsub(as.character(Rep_seq[1]), "", as.character(Rep_seq)) == ""){
        homo_tmp_2 = Pre_rep_status_tmp + length(Rep_seq)
      }
    }
  }
  Homopolymer_status = max(homo_tmp_1, homo_tmp_2)
  if((Post_rep_status + Post_rep_status) > 0 & gsub(as.character(Rep_B[1]), "", as.character(Rep_B)) == ""){
    Homopolymer_status = homo_tmp_1 + homo_tmp_2 - length(Rep_B)
  }
  return(list(Pre_rep_status, Post_rep_status, Homopolymer_status))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
