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
  short_support_length_adj <- NULL
  half_length <- NULL

  if (dim(msec)[1] > 0) {
    msec <- msec %>% mutate(
      distant_homology_rate = fun_zero(distant_homology, total_read)
    )
    msec <- msec %>% mutate(
      pre_minimum_length_adj = ifelse(indel_status == 1,
          (((pre_minimum_length - (indel_length + penalty_pre + 1)) +
          abs(pre_minimum_length - (indel_length + penalty_pre + 1))) / 2) +
            (indel_length + penalty_pre + 1),
          (((pre_minimum_length - penalty_pre) +
              abs(pre_minimum_length - penalty_pre)) / 2) +
            penalty_pre),
      post_minimum_length_adj = ifelse(indel_status == 1,
          (((post_minimum_length - (indel_length + penalty_post + 1)) +
          abs(post_minimum_length - (indel_length + penalty_post + 1))) / 2) +
            (indel_length + penalty_post + 1),
          (((post_minimum_length - penalty_post) +
              abs(post_minimum_length - penalty_post)) / 2) +
            penalty_post),
      pre_support_length_adj = ifelse(indel_status == 1,
          (((pre_support_length -
               (read_length - indel_length - penalty_post - 1)) -
        abs(pre_support_length -
              (read_length - indel_length - penalty_post - 1))) / 2) +
            (read_length - indel_length - penalty_post - 1),
        (((pre_support_length - (read_length - alt_length - penalty_post)) -
            abs(pre_support_length -
                  (read_length - alt_length - penalty_post))) / 2) +
          (read_length - alt_length - penalty_post)),
      post_support_length_adj = ifelse(indel_status == 1,
          (((post_support_length -
               (read_length - indel_length - penalty_pre - 1)) -
              abs(post_support_length -
                    (read_length - indel_length - penalty_pre - 1))) / 2) +
            (read_length - indel_length - penalty_pre - 1),
          (((post_support_length - (read_length - alt_length - penalty_pre)) -
              abs(post_support_length -
                    (read_length - alt_length - penalty_pre))) / 2) +
            (read_length - alt_length - penalty_pre))
    )
    msec <- msec %>% mutate(
      pre_minimum_length_adj = ifelse(indel_status == 1,
          (((pre_minimum_length_adj - pre_rep_status - 1) +
              abs(pre_minimum_length_adj - pre_rep_status - 1)) / 2) +
            pre_rep_status + 1,
          pre_minimum_length_adj),
      post_minimum_length_adj = ifelse(indel_status == 1,
          (((post_minimum_length_adj - post_rep_status - 1) +
              abs(post_minimum_length_adj - post_rep_status - 1)) / 2) +
            post_rep_status + 1,
          post_minimum_length_adj),
      pre_support_length_adj = ifelse(indel_status == 1,
       (((pre_support_length_adj - (read_length - post_rep_status - 1)) -
       abs(pre_support_length_adj - (read_length - post_rep_status - 1))) / 2) +
            (read_length - post_rep_status - 1),
          pre_support_length_adj),
      post_support_length_adj = ifelse(indel_status == 1,
       (((post_support_length_adj - (read_length - pre_rep_status - 1)) -
       abs(post_support_length_adj - (read_length - pre_rep_status - 1))) / 2) +
            (read_length - pre_rep_status - 1),
          post_support_length_adj)
    )
    msec <- msec %>% mutate(
      shortest_support_length_adj =
        (((pre_minimum_length_adj - post_minimum_length_adj) -
            abs(pre_minimum_length_adj - post_minimum_length_adj)) / 2) +
        post_minimum_length_adj,
      minimum_length_1 =
        (((pre_rep_status - (indel_length + penalty_pre + 1)) +
            abs(pre_rep_status - (indel_length + penalty_pre + 1))) / 2) +
        (indel_length + penalty_pre + 1),
      minimum_length_2 =
        (((post_rep_status - (indel_length + penalty_post + 1)) +
            abs(post_rep_status - (indel_length + penalty_post + 1))) / 2) +
        (indel_length + penalty_post + 1)
    )
    msec <- msec %>% mutate(
      minimum_length_1 = ifelse(indel_status == 1,
                                minimum_length_1, penalty_pre),
      minimum_length_2 = ifelse(indel_status == 1,
                                minimum_length_2, penalty_post)
    )
    msec <- msec %>% mutate(
      minimum_length =
        (((minimum_length_1 -  minimum_length_2) -
            abs(minimum_length_1 -  minimum_length_2)) / 2) +
        minimum_length_2
    )
    msec <- msec %>% mutate(
      shortest_support_length_adj =
        ((minimum_length - shortest_support_length_adj) +
        abs(minimum_length - shortest_support_length_adj)) / 2 +
        shortest_support_length_adj,
      pre_minimum_length_adj =
        ((minimum_length_1 - pre_minimum_length_adj) +
           abs(minimum_length_1 - pre_minimum_length_adj)) / 2 +
        pre_minimum_length_adj,
      post_minimum_length_adj =
        ((minimum_length_2 - post_minimum_length_adj) +
           abs(minimum_length_2 - post_minimum_length_adj)) / 2 +
        post_minimum_length_adj
    )
    msec <- msec %>% mutate(
      short_support_length_adj =
        (((short_support_length - minimum_length) +
            abs(short_support_length - minimum_length)) / 2) +
        minimum_length,
      altered_length = ifelse(mut_type == "snv",
                              alt_length,
                              ifelse(mut_type == "ins",
                                     indel_length,
                                     0)),
      half_length = as.integer((read_length - altered_length) / 2)
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
        read_length - altered_length - minimum_length_1 - minimum_length_2 + 1
    )
  }
  return(msec)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
