#' Summarizing function.
#'
#' This function summarizes the filtering results.
#'
#' @param msec Mutation filtering information.
#' @return msec
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @examples
#' data(msec_homology_searched)
#' fun_summary(msec_homology_searched)
#' @export
fun_summary <- function(msec) {
  distant_homology <- NULL
  total_read <- NULL
  indel_status <- NULL
  pre_minimum_length <- NULL
  indel_length <- NULL
  penalty_pre <- NULL
  post_minimum_length <- NULL
  penalty_post <- NULL
  pre_support_length <- NULL
  read_length <- NULL
  post_support_length <- NULL
  pre_minimum_length_adj <- NULL
  pre_rep_status <- NULL
  post_minimum_length_adj <- NULL
  post_rep_status <- NULL
  pre_support_length_adj <- NULL
  post_support_length_adj <- NULL
  minimum_length_1 <- NULL
  minimum_length_2 <- NULL
  minimum_length <- NULL
  shortest_support_length_adj <- NULL
  short_support_length <- NULL
  mut_type <- NULL
  alt_length <- NULL
  altered_length <- NULL
  altered_length2 <- NULL
  altered_length3 <- NULL
  short_support_length_adj <- NULL
  half_length <- NULL
  soft_clipped_rate <- NULL
  soft_clipped_read <- NULL
  
  if (dim(msec)[1] > 0) {
    msec <- msec %>% mutate(
      distant_homology_rate = fun_zero(distant_homology, total_read),
      altered_length = ifelse(mut_type == "snv",
                              alt_length,
                              ifelse(mut_type == "ins",
                                     indel_length,
                                     0)),
      altered_length2 = ifelse(mut_type == "del",
                              indel_length + 1,
                              alt_length),
      altered_length3 = ifelse(mut_type == "snv",
                               1,
                               ifelse(mut_type == "ins",
                                      1 + indel_length,
                                      1 - indel_length))
    )
    msec <- msec %>% mutate(
      penalty_pre = ifelse(pre_support_length == 0, 0, penalty_pre),
      penalty_post = ifelse(pre_support_length == 0, 0, penalty_post)
    )
    msec <- msec %>% mutate(
      pre_minimum_length_adj = (((pre_minimum_length - penalty_pre) +
              abs(pre_minimum_length - penalty_pre)) / 2) +
            penalty_pre,
      post_minimum_length_adj =
        (((post_minimum_length - penalty_post - altered_length2 + 1) +
          abs(post_minimum_length - penalty_post - altered_length2 + 1)) / 2) +
            penalty_post + altered_length2 - 1,
      pre_support_length_adj = 
        (((pre_support_length - (read_length - altered_length - penalty_post)) -
            abs(pre_support_length -
                  (read_length - altered_length - penalty_post))) / 2) +
          (read_length - altered_length - penalty_post),
      post_support_length_adj =
        (((post_support_length -
        (read_length - altered_length + altered_length2 - penalty_pre - 1)) -
            abs(post_support_length -
    (read_length - altered_length + altered_length2 - penalty_pre - 1))) / 2) +
          (read_length - altered_length + altered_length2 - penalty_pre - 1),
      short_support_length_adj =
        (short_support_length -
             (penalty_post + altered_length2 - 1) +
       abs(short_support_length -
             (penalty_post + altered_length2 - 1))) / 2 +
      (penalty_post + altered_length2 - 1)
    )
    msec <- msec %>% mutate(
      short_support_length_adj =
        (short_support_length_adj - penalty_pre +
           abs(short_support_length_adj - penalty_pre)) / 2 +
        penalty_pre
    )
    msec <- msec %>% mutate(
      shortest_support_length_adj =
        (((pre_minimum_length_adj - post_minimum_length_adj) -
            abs(pre_minimum_length_adj - post_minimum_length_adj)) / 2) +
        post_minimum_length_adj,
      minimum_length_1 = penalty_pre,
      minimum_length_2 = penalty_post + altered_length2 - 1
    )
    msec <- msec %>% mutate(
      minimum_length =
        (((minimum_length_1 -  minimum_length_2) -
            abs(minimum_length_1 -  minimum_length_2)) / 2) +
        minimum_length_2
    )
   msec <- msec %>% mutate(
      half_length = as.integer((read_length - altered_length3) / 2)
    )
    msec <- msec %>% mutate(
      short_support_length_total =
        short_support_length_adj - shortest_support_length_adj + 1,
      pre_support_length_total =
        pre_support_length_adj - pre_minimum_length_adj + 1,
      post_support_length_total =
        post_support_length_adj - post_minimum_length_adj + 1,
      half_length_total =
        half_length - minimum_length + 1,
      total_length_total =
        read_length - altered_length - minimum_length_1 - minimum_length_2 + 1,
      soft_clipped_rate = fun_zero(soft_clipped_read, total_read)
    )
  }
  return(msec)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
