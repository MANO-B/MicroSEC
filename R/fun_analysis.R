#' Analyzing function.
#'
#' This function analyzes the filtering results.
#'
#' @param MSEC Mutation filtering information.
#' @param Mut_depth Mutation coverage data.
#' @param Short_Homology_search_length Small sequence for homology search.
#' @param Minimum_Homology_search_length The sequence length for homology
#'   search.
#' @param threshold_p The largest p value of significant errors.
#' @param threshold_hairpin_ratio The smallest hairpin read ratio.
#' @param threshold_soft_clip_ratio The smallest rate of significantly 
#'   soft-clipped reads.
#' @param threshold_short_length Reads shorter than that are analyzed. 
#' @param threshold_distant_homology The smallest rate of reads from other 
#'   regions.
#' @param threshold_low_quality_rate The smallest rate of low quality bases.
#' @param Homopolymer_length The smallest length of homopolymers.
#' @return MSEC
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @export
fun_analysis = function(MSEC,
                        Mut_depth,
                        Short_Homology_search_length = 4,
                        Minimum_Homology_search_length = 40,
                        threshold_p = 10^(-6),
                        threshold_hairpin_ratio = 0.50,
                        threshold_soft_clip_ratio = 0.90,
                        threshold_short_length = 0.8,
                        threshold_distant_homology = 0.2,
                        threshold_low_quality_rate = 0.1,
                        Homopolymer_length = 15){
  if(!is.null(MSEC)){
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
               TRUE, FALSE),
      Not_long_repeat = 
        ifelse((Short_Homology_search_length +
                  Pre_rep_status +
                  Post_rep_status + 
                  Alt_length < Minimum_Homology_search_length),
               TRUE, FALSE)
    )
    MSEC$short_support_length_adjust_sum = 
      mapply(function(x, y) {return (Mut_depth[x,y])},
             1:dim(MSEC)[1], 
             MSEC$READ_length -
             MSEC$Altered_length + 2 -
             MSEC$shortest_support_length_adjust) -
      mapply(function(x, y) {return (Mut_depth[x,y])},
             1:dim(MSEC)[1], 
             MSEC$READ_length - 
             MSEC$Altered_length + 1 - 
             MSEC$short_support_length_adjust) + 
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1],
             MSEC$short_support_length_adjust + 2) -
      mapply(function(x, y) {return (Mut_depth[x,y])},
             1:dim(MSEC)[1],
             MSEC$shortest_support_length_adjust + 1)
    MSEC$Pre_support_length_adjust_sum =
      mapply(function(x, y) {return (Mut_depth[x,y])},
             1:dim(MSEC)[1],
             MSEC$Pre_support_length_adjust + 2) - 
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1], 
             MSEC$Pre_Minimum_length_adjust + 1)
    MSEC$Post_support_length_adjust_sum =
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1],
             MSEC$READ_length - 
             MSEC$Altered_length + 2 - 
             MSEC$Post_Minimum_length_adjust) - 
      mapply(function(x, y) {return (Mut_depth[x,y])},
             1:dim(MSEC)[1], 
             MSEC$READ_length - 
             MSEC$Altered_length - 
             MSEC$Post_support_length_adjust + 1)
    MSEC$Half_length_adjust_sum =
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1], 
             MSEC$READ_length - 
             MSEC$Altered_length + 2 - 
             MSEC$minimum_length) - 
      mapply(function(x, y) {return (Mut_depth[x,y])},
             1:dim(MSEC)[1],
             MSEC$READ_length - 
             MSEC$Altered_length - 
             MSEC$Half_length + 1) +
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1],
             MSEC$Half_length + 2) -
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1],
             MSEC$minimum_length + 1)
    MSEC$Total_length_adjust_sum = 
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1], 
             MSEC$READ_length - 
             MSEC$Altered_length + 2 - 
             MSEC$minimum_length_2) - 
      mapply(function(x, y) {return (Mut_depth[x,y])}, 
             1:dim(MSEC)[1], 
             MSEC$minimum_length_1 + 1)
    MSEC = MSEC %>% mutate(
      Short_short_support_sum = 
        (short_support_length_adjust_sum <= 
           Half_length_adjust_sum * threshold_short_length),
      Short_pre_support_sum = 
        (Pre_support_length_adjust_sum <= 
           Total_length_adjust_sum * threshold_short_length),
      Short_post_support_sum = 
        (Post_support_length_adjust_sum <= 
           Total_length_adjust_sum * threshold_short_length)
    )
    MSEC = MSEC %>% mutate(
      prob_Filter_1 = 
        fun_zero(short_support_length_adjust_sum, 
                 Half_length_adjust_sum) ^ Total_read,
      prob_Filter_3_pre = 
        fun_zero(Pre_support_length_adjust_sum, 
                 Total_length_adjust_sum) ^ Total_read,
      prob_Filter_3_post =
        fun_zero(Post_support_length_adjust_sum, 
                 Total_length_adjust_sum) ^ Total_read
    )
    MSEC = MSEC %>% mutate(
      prob_Filter_1 = ifelse((prob_Filter_1 > 1),
                             1, prob_Filter_1),
      prob_Filter_3_pre = ifelse((prob_Filter_3_pre > 1),
                                 1, prob_Filter_3_pre),
      prob_Filter_3_post = ifelse((prob_Filter_3_post > 1),
                                  1, prob_Filter_3_post)
    )
    MSEC = MSEC %>% mutate(
      Filter_1_mutation_intra_hairpin_loop = 
        ifelse((Short_short_support &
                  Short_short_support_sum &
                  prob_Filter_1 < threshold_p),
               TRUE, FALSE), 
      Filter_2_hairpin_structure = 
        ifelse((fun_zero(FLAG_Hairpin, Total_read) > threshold_hairpin_ratio),
               TRUE, FALSE),
      Filter_3_microhomology_induced_mutation = 
        ifelse((High_rate_Q18 & 
                  ((prob_Filter_3_pre <= threshold_p & 
                      Short_pre_support &
                      Short_pre_support_sum) | 
                   (prob_Filter_3_post <= threshold_p &
                      Short_post_support &
                      Short_post_support_sum))),
               TRUE, FALSE),
      Filter_4_soft_clipping = 
        ifelse((fun_zero(Soft_Clipped_read, Total_read)) >
                 threshold_soft_clip_ratio,
               TRUE, FALSE),
      Filter_5_highly_homologous_region =
        ifelse((distant_homology_rate >= threshold_distant_homology & 
                  Not_long_repeat),
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
      Caution = 
        ifelse((distant_homology_rate >= threshold_distant_homology & 
                  !Not_long_repeat),
               paste(Caution, 
                     "too repetitive to analyze homology,"),
               Caution)
    )
    MSEC = MSEC %>% mutate(
      Caution = 
        ifelse((prob_Filter_1 < threshold_p &
                !Filter_1_mutation_intra_hairpin_loop),
               paste(Caution,
                     "Filter 1: p is small, but supported enough long,"),
               Caution)
    )
    MSEC = MSEC %>% mutate(
      Caution = 
        ifelse(((prob_Filter_3_pre < threshold_p |
                  prob_Filter_3_post < threshold_p) &
                  !Filter_3_microhomology_induced_mutation),
               paste(Caution, 
                     "Filter 3: p is small, but supported enough long,"), 
               Caution)
    )
    MSEC = MSEC %>% mutate(
        MSEC_filter_1234 = ifelse(
          Filter_1_mutation_intra_hairpin_loop |
          Filter_2_hairpin_structure | 
          Filter_3_microhomology_induced_mutation | 
          Filter_4_soft_clipping,
          "Artifact suspicious", ""),
        MSEC_filter_12345 = ifelse(
          Filter_1_mutation_intra_hairpin_loop |
          Filter_2_hairpin_structure | 
          Filter_3_microhomology_induced_mutation | 
          Filter_4_soft_clipping |
          Filter_5_highly_homologous_region,
          "Artifact suspicious", ""),
        MSEC_filter_all = ifelse(
          Filter_1_mutation_intra_hairpin_loop |
          Filter_2_hairpin_structure | 
          Filter_3_microhomology_induced_mutation | 
          Filter_4_soft_clipping |
          Filter_5_highly_homologous_region | 
          Filter_6_simple_repeat |
          Filter_7_C_to_T_artifact |
          Filter_8_mutation_at_homopolymer,
          "Artifact suspicious", ""),
        Comment = Caution
      )
    MSEC = MSEC %>% select(-mut_type, -Alt_length, -Hairpin_length, 
            -Pre_Minimum_length, -Post_Minimum_length,
            -Low_quality_base_rate_under_Q18, -Pre_rep_status, -Post_rep_status, 
            -Homopolymer_status, -indel_status, -indel_length, -penalty_Pre,
            -penalty_Post, -Caution, -Pre_Minimum_length_adjust, -Half_length,
            -Post_Minimum_length_adjust, -Pre_support_length_adjust,
            -Post_support_length_adjust, -shortest_support_length_adjust,
            -minimum_length_1, -minimum_length_2, -minimum_length,
            -short_support_length_adjust, -Altered_length, -distant_homology,
            -short_support_length_total, -Pre_support_length_total,
            -Post_support_length_total, -Half_length_total,
            -Total_length_total, -Short_short_support, -Short_pre_support,
            -Short_post_support, -High_rate_Q18, -Not_long_repeat,
            -Short_short_support_sum, -Short_pre_support_sum,
            -Short_post_support_sum, -short_support_length_adjust_sum,
            -Pre_support_length_adjust_sum, -Post_support_length_adjust_sum,
            -Half_length_adjust_sum, -Total_length_adjust_sum
            )
  }
  return(MSEC)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
