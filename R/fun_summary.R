#' Summarizing function.
#'
#' This function summarizes the filtering results.
#'
#' @param msec Mutation filtering information.
#' @return msec
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @export
fun_summary <- function(msec) {
  if (!is.null(msec)) {
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
        (((pre_support_length - (read_length - penalty_post)) -
            abs(pre_support_length - (read_length - penalty_post))) / 2) +
          (read_length - penalty_post)),
      post_support_length_adj = ifelse(indel_status == 1,
          (((post_support_length -
               (read_length - indel_length - penalty_pre - 1)) -
              abs(post_support_length -
                    (read_length - indel_length - penalty_pre - 1))) / 2) +
            (read_length - indel_length - penalty_pre - 1),
          (((post_support_length - (read_length - penalty_pre)) -
              abs(post_support_length - (read_length - penalty_pre))) / 2) +
            (read_length - penalty_pre))
    )
    msec <- msec %>% mutate(
      pre_minimum_length_adj = ifelse(indel_status == 1,
          (((pre_minimum_length_adj - pre_rep_status) +
              abs(pre_minimum_length_adj - pre_rep_status)) / 2) +
            pre_rep_status,
          pre_minimum_length_adj),
      post_minimum_length_adj = ifelse(indel_status == 1,
          (((post_minimum_length_adj - post_rep_status) +
              abs(post_minimum_length_adj - post_rep_status)) / 2) +
            post_rep_status,
          post_minimum_length_adj),
      pre_support_length_adj = ifelse(indel_status == 1,
          (((pre_support_length_adj - (read_length - post_rep_status)) -
         abs(pre_support_length_adj - (read_length - post_rep_status))) / 2) +
            (read_length - post_rep_status),
          pre_support_length_adj),
      post_support_length_adj = ifelse(indel_status == 1,
          (((post_support_length_adj - (read_length - pre_rep_status)) -
         abs(post_support_length_adj - (read_length - pre_rep_status))) / 2) +
            (read_length - pre_rep_status),
          post_support_length_adj)
    )
    msec <- msec %>% mutate(
      shortest_support_length_adj =
        (((pre_minimum_length_adj - post_minimum_length_adj) -
            abs(pre_minimum_length_adj - post_minimum_length_adj)) / 2) +
        post_minimum_length_adj,
      minimum_length_1 =
        (((pre_rep_status - (indel_length + penalty_post + 1)) +
            abs(pre_rep_status - (indel_length + penalty_post + 1))) / 2) +
        (indel_length + penalty_post + 1),
      minimum_length_2 =
        (((post_rep_status - (indel_length + penalty_pre + 1)) +
            abs(post_rep_status - (indel_length + penalty_pre + 1))) / 2) +
        (indel_length + penalty_pre + 1)
    )
    msec <- msec %>% mutate(
      minimum_length_1 = ifelse(indel_status == 1, minimum_length_1, 0),
      minimum_length_2 = ifelse(indel_status == 1, minimum_length_2, 0)
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
