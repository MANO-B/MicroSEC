#' Summarizing function.
#'
#' This function summarizes the filtering results.
#'
#' @param MicroSEC Mutation filtering information.
#' @return MicroSEC
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom MicroSEC fun_zero
#' @export
fun_summary = function(MicroSEC){
  # thresholds
  threshold_p = 10^(-6)
  threshold_hairpin_ratio = 0.50
  threshold_hairpin_length = 30
  threshold_soft_clip_ratio = 0.90
  threshold_short_length = 0.8
  threshold_distant_homology = 0.2
  threshold_low_quality_rate = 0.1
  Homopolymer_length = 15
  
  MicroSEC = MicroSEC %>% mutate(
    distant_homology_rate = fun_zero(distant_homology, Total_read)
  )
  MicroSEC = MicroSEC %>% mutate(
    Pre_Minimum_length = ifelse(indel_status == 1, 
        (((Pre_Minimum_length - (indel_length + 1)) + 
            abs(Pre_Minimum_length - (indel_length + 1))) / 2) + 
          (indel_length + 1),
        Pre_Minimum_length),
    Post_Minimum_length = ifelse(indel_status == 1,
        (((Post_Minimum_length - (indel_length + 1)) +
            abs(Post_Minimum_length - (indel_length + 1))) / 2) +
          (indel_length + 1),
        Post_Minimum_length),
    Pre_support_length = ifelse(indel_status == 1, 
        (((Pre_support_length - (READ_length - indel_length - 1)) - 
            abs(Pre_support_length - (READ_length - indel_length - 1))) / 2) +
          (READ_length - indel_length - 1), 
        Pre_support_length),
    Post_support_length = ifelse(indel_status == 1,
        (((Post_support_length - (READ_length - indel_length - 1)) -
            abs(Post_support_length - (READ_length - indel_length - 1))) / 2) +
          (READ_length - indel_length - 1),
        Post_support_length)
  )
  MicroSEC = MicroSEC %>% mutate(
    Pre_Minimum_length = ifelse(indel_status == 1, 
        (((Pre_Minimum_length - Pre_rep_status) + 
            abs(Pre_Minimum_length - Pre_rep_status)) / 2) + 
          Pre_rep_status, 
        Pre_Minimum_length),
    Post_Minimum_length = ifelse(indel_status == 1,
        (((Post_Minimum_length - Post_rep_status) + 
            abs(Post_Minimum_length -Post_rep_status)) / 2) + 
          Post_rep_status,
        Post_Minimum_length),
    Pre_support_length = ifelse(indel_status == 1,
        (((Pre_support_length - (READ_length - Post_rep_status)) -
            abs(Pre_support_length - (READ_length - Post_rep_status))) / 2) +
          (READ_length - Post_rep_status),
        Pre_support_length),
    Post_support_length = ifelse(indel_status == 1, 
        (((Post_support_length - (READ_length - Pre_rep_status)) -
            abs(Post_support_length - (READ_length - Pre_rep_status))) / 2) +
          (READ_length - Pre_rep_status),
        Post_support_length)
  )
  MicroSEC = MicroSEC %>% mutate(
    shortest_support_length =
      (((Pre_Minimum_length - Post_Minimum_length) - 
          abs(Pre_Minimum_length - Post_Minimum_length)) / 2) + 
      Post_Minimum_length,
    minimum_length_1 = 
      (((Pre_rep_status - (indel_length + 1)) + 
          abs(Pre_rep_status - (indel_length + 1))) / 2) + 
      (indel_length + 1),
    minimum_length_2 =
      (((Post_rep_status - (indel_length + 1)) + 
          abs(Post_rep_status - (indel_length + 1))) / 2) + 
      (indel_length + 1)
  )
  MicroSEC = MicroSEC %>% mutate(
    minimum_length_1 = ifelse(indel_status == 1, minimum_length_1, 0),
    minimum_length_2 = ifelse(indel_status == 1, minimum_length_2, 0)
  )
  MicroSEC = MicroSEC %>% mutate(
    minimum_length =
      (((minimum_length_1 -  minimum_length_2) - 
          abs(minimum_length_1 -  minimum_length_2)) / 2) + 
      minimum_length_2
  )
  MicroSEC = MicroSEC %>% mutate(
    short_support_length =
      (((short_support_length - minimum_length) + 
          abs(short_support_length - minimum_length)) / 2) + 
      minimum_length,
    Altered_length = ifelse(mut_type == "snv", 
                            Alt_length,
                            ifelse(mut_type == "ins",
                                   indel_length,
                                   0)),
  )
  MicroSEC = MicroSEC %>% mutate(
    short_support_length_adjust = 
      short_support_length - shortest_support_length + 1,
    Pre_support_length_adjust =
      Pre_support_length - Pre_Minimum_length + 1,
    Post_support_length_adjust = 
      Post_support_length - Post_Minimum_length + 1,
    Half_length_adjust =
      as.integer((READ_length - Altered_length) / 2)  - minimum_length + 1,
    Total_length_adjust = 
      READ_length - Altered_length - minimum_length_1 - minimum_length_2 + 1
  )
  MicroSEC = MicroSEC %>% select(
    -minimum_length_1, 
    -minimum_length_2,
    -minimum_length, 
    -Altered_length, 
    -shortest_support_length)
  MicroSEC = MicroSEC %>% mutate(
    Short_short_support = 
      (short_support_length_adjust <= 
         threshold_short_length * Half_length_adjust),
    Short_pre_support = 
      (Pre_support_length_adjust <= 
         threshold_short_length * Total_length_adjust),
    Short_post_support = 
      (Post_support_length_adjust <= 
         threshold_short_length * Total_length_adjust),
    High_rate_Q18 = 
      ifelse((Low_quality_base_rate_under_Q18 < threshold_low_quality_rate),
             TRUE, FALSE),
    prob_Filter_1 = 
      (short_support_length_adjust / Half_length_adjust) ^ Total_read,
    prob_Filter_3_pre = 
      (Pre_support_length_adjust / Total_length_adjust) ^ Total_read,
    prob_Filter_3_post =
      (Post_support_length_adjust / Total_length_adjust) ^ Total_read
  )
  MicroSEC = MicroSEC %>% mutate(
    prob_Filter_1 = ifelse((prob_Filter_1 > 1), 1, prob_Filter_1),
    prob_Filter_3_pre = ifelse((prob_Filter_3_pre > 1), 1, prob_Filter_3_pre),
    prob_Filter_3_post = ifelse((prob_Filter_3_post > 1), 1, prob_Filter_3_post)
  )
  MicroSEC = MicroSEC %>% mutate(
    Filter_1_mutation_intra_hairpin_loop = 
      ifelse((Short_short_support & prob_Filter_1 < threshold_p),
             TRUE, FALSE), 
    Filter_2_hairpin_structure = 
      ifelse((zero(FLAG_Hairpin, Total_read) > threshold_hairpin_ratio |
                Hairpin_length >= threshold_hairpin_length),
             TRUE, FALSE),
    Filter_3_microhomology_induced_mutation = 
      ifelse((High_rate_Q18 & 
                ((prob_Filter_3_pre <= threshold_p & Short_pre_support) | 
                   (prob_Filter_3_post <= threshold_p & Short_post_support))),
             TRUE, FALSE),
    Filter_4_soft_clipping = 
      ifelse((zero(Soft_Clipped_read, Total_read)) > threshold_soft_clip_ratio,
             TRUE, FALSE),
    Filter_5_highly_homologous_region =
      ifelse((distant_homology_rate >= threshold_distant_homology),
             TRUE, FALSE),
    Filter_6_simple_repeat =
      ifelse((SimpleRepeat_TRF == "Y"),
             TRUE, FALSE), 
    Filter_7_C_to_T_artifact =
      ifelse((Transition == "C>T_g_FFPE"), 
             TRUE, FALSE),
    Filter_8_mutation_at_homopolymer = 
      ifelse((Homopolymer_status  >= Homopolymer_length), 
             TRUE, FALSE)
  )
  MicroSEC = MicroSEC %>% mutate(
    MicroSEC_filter_1234 = 
      Filter_1_mutation_intra_hairpin_loop |
      Filter_2_hairpin_structure | 
      Filter_3_microhomology_induced_mutation | 
      Filter_4_soft_clipping,
    MicroSEC_filter_12345 = 
      MicroSEC_filter_1234 | 
      Filter_5_highly_homologous_region,
    MicroSEC_filter_all = 
      MicroSEC_filter_12345 | 
      Filter_6_simple_repeat |
      Filter_7_C_to_T_artifact |
      Filter_8_mutation_at_homopolymer
  )
  return(MicroSEC)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
