#' Summarizing function.
#'
#' This function summarizes the filtering results.
#'
#' @param MSEC Mutation filtering information.
#' @return MSEC
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @export
fun_summary = function(MSEC){
  MSEC = MSEC %>% mutate(
    distant_homology_rate = fun_zero(distant_homology, Total_read)
  )
  MSEC = MSEC %>% mutate(
    Pre_Minimum_length_adjust = ifelse(indel_status == 1, 
        (((Pre_Minimum_length - (indel_length + penalty_Pre + 1)) + 
            abs(Pre_Minimum_length - (indel_length + penalty_Pre + 1))) / 2) + 
          (indel_length + penalty_Pre + 1),
        (((Pre_Minimum_length - penalty_Pre) + 
            abs(Pre_Minimum_length - penalty_Pre)) / 2) + 
          penalty_Pre),
    Post_Minimum_length_adjust = ifelse(indel_status == 1,
        (((Post_Minimum_length - (indel_length + penalty_Post + 1)) +
            abs(Post_Minimum_length - (indel_length + penalty_Post + 1))) / 2) +
          (indel_length + penalty_Post + 1),
        (((Post_Minimum_length - penalty_Post) + 
            abs(Post_Minimum_length - penalty_Post)) / 2) + 
          penalty_Post),
    Pre_support_length_adjust = ifelse(indel_status == 1, 
        (((Pre_support_length - (READ_length - indel_length - penalty_Post - 1)) - 
      abs(Pre_support_length - (READ_length - indel_length - penalty_Post - 1))) / 2) +
          (READ_length - indel_length - penalty_Post - 1), 
      (((Pre_support_length - (READ_length - penalty_Post)) - 
          abs(Pre_support_length - (READ_length - penalty_Post))) / 2) +
        (READ_length - penalty_Post)),
    Post_support_length_adjust = ifelse(indel_status == 1,
        (((Post_support_length - (READ_length - indel_length - penalty_Pre - 1)) -
            abs(Post_support_length - (READ_length - indel_length - penalty_Pre - 1))) / 2) +
          (READ_length - indel_length - penalty_Pre - 1),
        (((Post_support_length - (READ_length - penalty_Pre)) -
            abs(Post_support_length - (READ_length - penalty_Pre))) / 2) +
          (READ_length - penalty_Pre))
  )
  MSEC = MSEC %>% mutate(
    Pre_Minimum_length_adjust = ifelse(indel_status == 1, 
        (((Pre_Minimum_length_adjust - Pre_rep_status) + 
            abs(Pre_Minimum_length_adjust - Pre_rep_status)) / 2) + 
          Pre_rep_status, 
        Pre_Minimum_length_adjust),
    Post_Minimum_length_adjust = ifelse(indel_status == 1,
        (((Post_Minimum_length_adjust - Post_rep_status) + 
            abs(Post_Minimum_length_adjust -Post_rep_status)) / 2) + 
          Post_rep_status,
        Post_Minimum_length_adjust),
    Pre_support_length_adjust = ifelse(indel_status == 1,
        (((Pre_support_length_adjust - (READ_length - Post_rep_status)) -
       abs(Pre_support_length_adjust - (READ_length - Post_rep_status))) / 2) +
          (READ_length - Post_rep_status),
        Pre_support_length_adjust),
    Post_support_length_adjust = ifelse(indel_status == 1, 
        (((Post_support_length_adjust - (READ_length - Pre_rep_status)) -
       abs(Post_support_length_adjust - (READ_length - Pre_rep_status))) / 2) +
          (READ_length - Pre_rep_status),
        Post_support_length_adjust)
  )
  MSEC = MSEC %>% mutate(
    shortest_support_length_adjust =
      (((Pre_Minimum_length_adjust - Post_Minimum_length_adjust) - 
          abs(Pre_Minimum_length_adjust - Post_Minimum_length_adjust)) / 2) + 
      Post_Minimum_length_adjust,
    minimum_length_1 = 
      (((Pre_rep_status - (indel_length + penalty_Post + 1)) + 
          abs(Pre_rep_status - (indel_length + penalty_Post + 1))) / 2) + 
      (indel_length + penalty_Post + 1),
    minimum_length_2 =
      (((Post_rep_status - (indel_length + penalty_Pre + 1)) + 
          abs(Post_rep_status - (indel_length + penalty_Pre + 1))) / 2) + 
      (indel_length + penalty_Pre + 1)
  )
  MSEC = MSEC %>% mutate(
    minimum_length_1 = ifelse(indel_status == 1, minimum_length_1, 0),
    minimum_length_2 = ifelse(indel_status == 1, minimum_length_2, 0)
  )
  MSEC = MSEC %>% mutate(
    minimum_length =
      (((minimum_length_1 -  minimum_length_2) - 
          abs(minimum_length_1 -  minimum_length_2)) / 2) + 
      minimum_length_2
  )
  MSEC = MSEC %>% mutate(
    short_support_length_adjust =
      (((short_support_length - minimum_length) + 
          abs(short_support_length - minimum_length)) / 2) + 
      minimum_length,
    Altered_length = ifelse(mut_type == "snv", 
                            Alt_length,
                            ifelse(mut_type == "ins",
                                   indel_length,
                                   0)),
    Half_length = as.integer((READ_length - Altered_length) / 2)
  )
  MSEC = MSEC %>% mutate(
    short_support_length_total = 
      short_support_length_adjust - shortest_support_length + 1,
    Pre_support_length_total =
      Pre_support_length_adjust - Pre_Minimum_length_adjust + 1,
    Post_support_length_total = 
      Post_support_length_adjust - Post_Minimum_length_adjust + 1,
    Half_length_total =
      Half_length - minimum_length + 1,
    Total_length_total = 
      READ_length - Altered_length - minimum_length_1 - minimum_length_2 + 1
  )
  return(MSEC)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
