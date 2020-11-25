#' Repeat check function.
#'
#' This function attempts to check the repetitive sequence around the mutation.
#'
#' @param rep_a The shorter sequence of Ref and Alt.
#' @param rep_b The longer sequence of Ref and Alt.
#' @param ref_seq Reference sequence around the mutation.
#' @param ref_width Search length for ref_seq.
#' @param del Insertion: 0, Deletion: 1
#' @return list(pre_rep_status, post_rep_status, homopolymer_status)
#' @importFrom Biostrings DNAString
#' @importFrom BiocGenerics as.data.frame
fun_repeat_check <- function(rep_a, rep_b, ref_seq, ref_width, del) {
  pre_rep_status <- 0
  post_rep_status <- 0
  homopolymer_status <- 0
  rep_b <- rep_b[2:length(rep_b)]
  homo_tmp_1 <- 0
  homo_tmp_2 <- 0
  for (q in seq_len(length(rep_b))) {
    rep_seq <- rep_b[q:length(rep_b)]
    check_rep <- TRUE
    post_rep_status_tmp <- 0
    for (r in (1 + del * length(rep_b)):(ref_width - 1)) {
      if (check_rep &
          ref_seq[ref_width + r + 1] ==
          rep_seq[((r - 1) %% length(rep_seq) + 1)]) {
        post_rep_status_tmp <- post_rep_status_tmp + 1
      }
      else {
        check_rep <- FALSE
      }
    }
    if (post_rep_status_tmp >= length(rep_seq)) {
      post_rep_status <- max(post_rep_status, post_rep_status_tmp)
      if (gsub(as.character(rep_seq[1]), "", as.character(rep_seq)) == "") {
        homo_tmp_1 <- post_rep_status_tmp + length(rep_seq)
      }
    }
  }
  for (q in seq_len(length(rep_b))) {
    rep_seq <- rep_b[1:q]
    check_rep <- TRUE
    pre_rep_status_tmp <- 0
    for (r in 1:(ref_width + 1)) {
      if (check_rep &
          ref_seq[ref_width + 2 - r] ==
            rep_seq[((length(rep_seq) - r) %% length(rep_seq) + 1)]) {
        pre_rep_status_tmp <- pre_rep_status_tmp + 1
      }
      else {
        check_rep <- FALSE
      }
    }
    if (pre_rep_status_tmp >= length(rep_seq)) {
      pre_rep_status <- max(pre_rep_status, pre_rep_status_tmp)
      if (gsub(as.character(rep_seq[1]), "", as.character(rep_seq)) == "") {
        homo_tmp_2 <- pre_rep_status_tmp + length(rep_seq)
      }
    }
  }
  homopolymer_status <- max(homo_tmp_1, homo_tmp_2)
  if ((post_rep_status + post_rep_status) > 0 &
      gsub(as.character(rep_b[1]), "", as.character(rep_b)) == "") {
    homopolymer_status <- homo_tmp_1 + homo_tmp_2 - length(rep_b)
  }
  return(list(pre_rep_status, post_rep_status, homopolymer_status))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
