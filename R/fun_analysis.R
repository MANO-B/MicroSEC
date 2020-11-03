#' Analyzing function.
#'
#' This function analyzes the filtering results.
#'
#' @param MSEC Mutation filtering information.
#' @return MSEC
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @export
fun_analysis = function(MSEC,
                        Mut_depth,
                        threshold_p = 10^(-6),
                        threshold_hairpin_ratio = 0.50,
                        threshold_hairpin_length = 30,
                        threshold_soft_clip_ratio = 0.90,
                        threshold_short_length = 0.8,
                        threshold_distant_homology = 0.2,
                        threshold_low_quality_rate = 0.1,
                        Homopolymer_length = 15){
  MSEC = MSEC %>% mutate(
    Short_short_support = 
      (short_support_length_total <= 
         threshold_short_length * Half_length_total),
    Short_pre_support = 
      (Pre_support_length_total <= 
         threshold_short_length * Total_length_total),
    Short_post_support = 
      (Post_support_length_total <= 
         threshold_short_length * Total_length_total),
    High_rate_Q18 = 
      ifelse((Low_quality_base_rate_under_Q18 < threshold_low_quality_rate),
             TRUE, FALSE)
  )
  MSEC$short_support_length_adjust_sum = 
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$READ_length + 1 - MSEC$short_support_length_adjust) - 
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$short_support_length_adjust + 1) -
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$READ_length + 1 - MSEC$shortest_support_length_adjust) +
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$shortest_support_length_adjust + 1)
  MSEC$Pre_support_length_adjust_sum =
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$Pre_support_length_adjust + 1) - 
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$Pre_Minimum_length_adjust + 1)
  MSEC$Post_support_length_adjust_sum =
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$READ_length + 1 - MSEC$Post_Minimum_length_adjust) - 
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$READ_length + 1 - MSEC$Post_support_length_adjust)
  MSEC$Half_length_adjust_sum =
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$READ_length + 1 - MSEC$minimum_length) - 
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$READ_length + 1 - MSEC$Half_length) -
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$Half_length + 1) +
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$minimum_length + 1)
  MSEC$Total_length_adjust_sum = 
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$READ_length + 1 - MSEC$minimum_length_2) - 
    mapply(function(x, y) {return (Mut_depth[x,y])}, 1:dim(MSEC)[1], MSEC$minimum_length_1 + 1)
  MSEC = MSEC %>% mutate(
    prob_Filter_1 = 
      (short_support_length_adjust_sum / Half_length_adjust_sum) ^ Total_read,
    prob_Filter_3_pre = 
      (Pre_support_length_adjust_sum / Total_length_adjust_sum) ^ Total_read,
    prob_Filter_3_post =
      (Post_support_length_adjust_sum / Total_length_adjust_sum) ^ Total_read
  )
  MSEC = MSEC %>% mutate(
    prob_Filter_1 = ifelse((prob_Filter_1 > 1), 1, prob_Filter_1),
    prob_Filter_3_pre = ifelse((prob_Filter_3_pre > 1), 1, prob_Filter_3_pre),
    prob_Filter_3_post = ifelse((prob_Filter_3_post > 1), 1, prob_Filter_3_post)
  )
  MSEC = MSEC %>% mutate(
    Filter_1_mutation_intra_hairpin_loop = 
      ifelse((Short_short_support & prob_Filter_1 < threshold_p),
             TRUE, FALSE), 
    Filter_2_hairpin_structure = 
      ifelse((fun_zero(FLAG_Hairpin, Total_read) > threshold_hairpin_ratio |
                Hairpin_length >= threshold_hairpin_length),
             TRUE, FALSE),
    Filter_3_microhomology_induced_mutation = 
      ifelse((High_rate_Q18 & 
                ((prob_Filter_3_pre <= threshold_p & Short_pre_support) | 
                   (prob_Filter_3_post <= threshold_p & Short_post_support))),
             TRUE, FALSE),
    Filter_4_soft_clipping = 
      ifelse((fun_zero(Soft_Clipped_read, Total_read)) > threshold_soft_clip_ratio,
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
  MSEC = MSEC %>% mutate(
    MSEC_filter_1234 = 
      Filter_1_mutation_intra_hairpin_loop |
      Filter_2_hairpin_structure | 
      Filter_3_microhomology_induced_mutation | 
      Filter_4_soft_clipping,
    MSEC_filter_12345 = 
      MSEC_filter_1234 | 
      Filter_5_highly_homologous_region,
    MSEC_filter_all = 
      MSEC_filter_12345 | 
      Filter_6_simple_repeat |
      Filter_7_C_to_T_artifact |
      Filter_8_mutation_at_homopolymer
  )
  return(MSEC)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
