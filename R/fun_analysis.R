#' Analyzing function.
#'
#' This function analyzes the filtering results.
#'
#' @param msec Mutation filtering information.
#' @param mut_depth Mutation coverage data.
#' @param short_homology_search_length Small sequence for homology search.
#' @param min_homology_search The sequence length for homology
#'   search.
#' @param threshold_p The largest p value of significant errors.
#' @param threshold_hairpin_ratio The smallest hairpin read ratio.
#' @param threshold_short_length Reads shorter than that are analyzed.
#' @param threshold_distant_homology The smallest rate of reads from other
#'   regions.
#' @param threshold_soft_clip_ratio The rate of soft-clipped reads.
#' @param threshold_low_quality_rate The smallest rate of low quality bases.
#' @param homopolymer_length The smallest length of homopolymers.
#' @return msec
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom BiocGenerics mapply
#' @examples
#' data(msec_summarized)
#' data(mut_depth_checked)
#' fun_analysis(msec = msec_summarized,
#'              mut_depth = mut_depth_checked,
#'              short_homology_search_length = 4,
#'              min_homology_search = 40,
#'              threshold_p = 10 ^ (-6),
#'              threshold_hairpin_ratio = 0.50,
#'              threshold_short_length = 0.75,
#'              threshold_distant_homology = 0.15,
#'              threshold_soft_clip_ratio = 0.50,
#'              threshold_low_quality_rate = 0.1,
#'              homopolymer_length = 15
#' )
#' @export
fun_analysis <- function(msec,
                        mut_depth,
                        short_homology_search_length,
                        min_homology_search ,
                        threshold_p,
                        threshold_hairpin_ratio,
                        threshold_short_length,
                        threshold_distant_homology,
                        threshold_soft_clip_ratio,
                        threshold_low_quality_rate,
                        homopolymer_length) {
  if (dim(msec)[1] > 0) {
    short_support_length_total <- NULL
    half_length_total <- NULL
    pre_support_length_total <- NULL
    total_length_total <- NULL
    post_support_length_total <- NULL
    low_quality_base_rate_under_q18 <- NULL
    pre_rep_status <- NULL
    post_rep_status <- NULL
    alt_length <- NULL
    short_support_length_adj_sum <- NULL
    half_length_adj_sum <- NULL
    pre_support_length_adj_sum <- NULL
    total_length_pre_adj_sum <- NULL
    total_length_post_adj_sum <- NULL
    post_support_length_adj_sum <- NULL
    total_read <- NULL
    prob_filter_1 <- NULL
    prob_filter_3_pre <- NULL
    prob_filter_3_post <- NULL
    short_short_support <- NULL
    short_short_support_sum <- NULL
    flag_hairpin <- NULL
    high_rate_q18 <- NULL
    short_pre_support <- NULL
    short_pre_support_sum <- NULL
    short_post_support <- NULL
    short_post_support_sum <- NULL
    soft_clipped_read <- NULL
    pre_farthest <- NULL
    post_farthest <- NULL    
    distant_homology_rate <- NULL
    not_long_repeat <- NULL
    SimpleRepeat_TRF <- NULL
    homopolymer_status <- NULL
    soft_clipped_rate <- NULL
    caution <- NULL
    filter_1_mutation_intra_hairpin_loop <- NULL
    filter_2_hairpin_structure <- NULL
    filter_3_microhomology_induced_mutation <- NULL
    filter_4_highly_homologous_region <- NULL
    filter_5_soft_clipped_reads <- NULL
    filter_6_simple_repeat <- NULL
    filter_7_mutation_at_homopolymer <- NULL
    filter_8_low_quality <- NULL
    mut_type <- NULL
    hairpin_length <- NULL
    pre_minimum_length <- NULL
    post_minimum_length<- NULL
    indel_status <- NULL
    indel_length <- NULL
    penalty_pre <- NULL
    penalty_post <- NULL
    pre_minimum_length_adj <- NULL
    half_length <- NULL
    post_minimum_length_adj <- NULL
    pre_support_length_adj <- NULL
    post_support_length_adj <- NULL
    shortest_support_length_adj <- NULL
    minimum_length_1 <- NULL
    minimum_length_2 <- NULL
    minimum_length <- NULL
    short_support_length_adj <- NULL
    altered_length <- NULL
    distant_homology <- NULL
    low_quality_pre <- NULL
    low_quality_post <- NULL
    
    mut_depth_pre <- mut_depth[[1]]
    mut_depth_post <- mut_depth[[2]]
    mut_depth_short <- mut_depth[[3]]
    
    
    msec <- msec %>% mutate(
      short_short_support =
        (short_support_length_total <=
           threshold_short_length * half_length_total),
      short_pre_support =
        (pre_support_length_total <=
           threshold_short_length * total_length_total),
      short_post_support =
        (post_support_length_total <=
           threshold_short_length * total_length_total),
      high_rate_q18 =
        ifelse((low_quality_base_rate_under_q18 < threshold_low_quality_rate),
               TRUE, FALSE),
      not_long_repeat =
        ifelse((short_homology_search_length +
                  pre_rep_status +
                  post_rep_status +
                  alt_length < min_homology_search),
               TRUE, FALSE)
    )
    msec$short_support_length_adj_sum <-
      mapply(
        function(x, y) {
          return(mut_depth_short[x, y])
        },
        1:dim(msec)[1],
        msec$short_support_length_adj + 2
      ) -
      mapply(
        function(x, y) {
          return(mut_depth_short[x, y])
        },
        1:dim(msec)[1],
        msec$shortest_support_length_adj + 2
      )
    msec$pre_support_length_adj_sum <-
      mapply(
        function(x, y) {
          return(mut_depth_pre[x, y])
        },
        1:dim(msec)[1],
        msec$pre_support_length_adj + 2
      ) -
      mapply(
        function(x, y) {
          return(mut_depth_pre[x, y])
        },
        1:dim(msec)[1],
        msec$pre_minimum_length_adj + 2
      )
    msec$post_support_length_adj_sum <-
      mapply(
        function(x, y) {
          return(mut_depth_post[x, y])
        },
        1:dim(msec)[1],
        msec$post_support_length_adj + 2
      ) -
      mapply(
        function(x, y) {
          return(mut_depth_post[x, y])
        },
        1:dim(msec)[1],
        msec$post_minimum_length_adj + 2
      )
    msec$half_length_adj_sum <-
      mapply(
        function(x, y) {
          return(mut_depth_short[x, y])
        },
        1:dim(msec)[1],
        msec$half_length + 2
      ) -
      mapply(
        function(x, y) {
          return(mut_depth_short[x, y])
        },
        1:dim(msec)[1],
        msec$minimum_length + 2
      )
    msec$total_length_pre_adj_sum <-
      mapply(
        function(x, y) {
          return(mut_depth_pre[x, y])
        },
        1:dim(msec)[1],
        msec$read_length -
          msec$altered_length + 2 -
          msec$minimum_length_2
      ) -
      mapply(
        function(x, y) {
          return(mut_depth_pre[x, y])
        },
        1:dim(msec)[1],
        msec$minimum_length_1 + 2
      )
    msec$total_length_post_adj_sum <-
      mapply(
        function(x, y) {
          return(mut_depth_post[x, y])
        },
        1:dim(msec)[1],
        msec$read_length -
          msec$altered_length + 2 -
          msec$minimum_length_2
      ) -
      mapply(
        function(x, y) {
          return(mut_depth_post[x, y])
        },
        1:dim(msec)[1],
        msec$minimum_length_1 + 2
      )
    msec <- msec %>% mutate(
      short_short_support_sum =
        (short_support_length_adj_sum <=
           half_length_adj_sum * threshold_short_length),
      short_pre_support_sum =
        (pre_support_length_adj_sum <=
           total_length_pre_adj_sum * threshold_short_length),
      short_post_support_sum =
        (post_support_length_adj_sum <=
           total_length_post_adj_sum * threshold_short_length)
    )
    msec <- msec %>% mutate(
      prob_filter_1 =
        fun_zero(short_support_length_adj_sum,
                 half_length_adj_sum) ^ total_read,
      prob_filter_3_pre =
        fun_zero(pre_support_length_adj_sum,
                 total_length_pre_adj_sum) ^ total_read,
      prob_filter_3_post =
        fun_zero(post_support_length_adj_sum,
                 total_length_post_adj_sum) ^ total_read
    )
    msec <- msec %>% mutate(
      prob_filter_1 = ifelse((prob_filter_1 > 1),
                             1, prob_filter_1),
      prob_filter_3_pre = ifelse((prob_filter_3_pre > 1),
                                 1, prob_filter_3_pre),
      prob_filter_3_post = ifelse((prob_filter_3_post > 1),
                                  1, prob_filter_3_post)
    )
    msec <- msec %>% mutate(
      filter_1_mutation_intra_hairpin_loop =
        ifelse((short_short_support &
                  short_short_support_sum &
                  prob_filter_1 < threshold_p),
               TRUE, FALSE),
      filter_2_hairpin_structure =
        ifelse((fun_zero(flag_hairpin, total_read) > threshold_hairpin_ratio),
               TRUE, FALSE),
      filter_3_microhomology_induced_mutation =
        ifelse((high_rate_q18 &
                  ((prob_filter_3_pre <= threshold_p &
                      short_pre_support &
                      short_pre_support_sum) |
                   (prob_filter_3_post <= threshold_p &
                      short_post_support &
                      short_post_support_sum))),
               TRUE, FALSE),
      filter_4_highly_homologous_region =
        ifelse((distant_homology_rate >= threshold_distant_homology &
                #(pre_farthest == pre_support_length | post_farthest == post_support_length) &
                  not_long_repeat),
               TRUE, FALSE),
      filter_5_soft_clipped_reads =
        ifelse((soft_clipped_rate >= threshold_soft_clip_ratio),
               TRUE, FALSE),
      filter_6_simple_repeat =
        ifelse((SimpleRepeat_TRF == "Y"),
               TRUE, FALSE),
      filter_7_mutation_at_homopolymer =
        ifelse((homopolymer_status  >= homopolymer_length),
               TRUE, FALSE),
      filter_8_low_quality =
        ifelse((low_quality_base_rate_under_q18 >= threshold_low_quality_rate |
                low_quality_pre >= threshold_low_quality_rate |
                low_quality_post >= threshold_low_quality_rate),
               TRUE, FALSE)
    )
    msec <- msec %>% mutate(
      caution =
        ifelse((distant_homology_rate >= threshold_distant_homology &
                  !not_long_repeat),
               paste(caution,
                     "too repetitive to analyze homology,"),
               caution)
    )
    msec <- msec %>% mutate(
      caution =
        ifelse((prob_filter_1 < threshold_p &
                  (!short_short_support |
                     !short_short_support_sum)),
               paste(caution,
                     "filter 1: p is small, but supported enough long,"),
               caution)
    )
    msec <- msec %>% mutate(
      caution =
        ifelse(((prob_filter_3_pre < threshold_p |
                  prob_filter_3_post < threshold_p) &
                  !filter_3_microhomology_induced_mutation &
                  high_rate_q18),
               paste(caution,
                     "filter 3: p is small, but supported enough long,"),
               caution)
    )
    msec <- msec %>% mutate(
      caution =
        ifelse(((prob_filter_3_pre < threshold_p |
                   prob_filter_3_post < threshold_p) &
                  !filter_3_microhomology_induced_mutation &
                  !high_rate_q18),
               paste(caution,
                     "filter 3: p is small, but reads are low quality,"),
               caution)
    )
    msec <- msec %>% mutate(
      caution =
        ifelse((distant_homology_rate >= threshold_distant_homology &
                  !filter_4_highly_homologous_region &
                  !not_long_repeat),
               paste(caution,
                     "filter 4: sequence is too repetitive,"),
               caution)
    )
    msec <- msec %>% mutate(
        msec_filter_123 = ifelse(
          filter_1_mutation_intra_hairpin_loop |
          filter_2_hairpin_structure |
          filter_3_microhomology_induced_mutation,
          "Artifact suspicious", ""),
        msec_filter_1234 = ifelse(
          filter_1_mutation_intra_hairpin_loop |
          filter_2_hairpin_structure |
          filter_3_microhomology_induced_mutation |
          filter_4_highly_homologous_region,
          "Artifact suspicious", ""),
        msec_filter_all = ifelse(
          filter_1_mutation_intra_hairpin_loop |
          filter_2_hairpin_structure |
          filter_3_microhomology_induced_mutation |
          filter_4_highly_homologous_region |
          filter_5_soft_clipped_reads |
          filter_6_simple_repeat |
          filter_7_mutation_at_homopolymer |
          filter_8_low_quality,
          "Artifact suspicious", ""),
        comment = caution
      )
    msec <- msec %>% select(-caution)
  }
  msec <- msec %>% select(-mut_type, -alt_length, -hairpin_length,
                          -pre_minimum_length, -post_minimum_length,
                          -pre_rep_status, -post_rep_status,
                          -homopolymer_status, -indel_status, -indel_length,
                          -penalty_pre,
                          -penalty_post, -caution, -pre_minimum_length_adj,
                          -half_length,
                          -post_minimum_length_adj, -pre_support_length_adj,
                          -post_support_length_adj,
                          -shortest_support_length_adj,
                          -minimum_length_1, -minimum_length_2, -minimum_length,
                          -short_support_length_adj, -altered_length,
                          -distant_homology,
                          -short_support_length_total,
                          -pre_support_length_total,
                          -post_support_length_total, -half_length_total,
                          -total_length_total, -high_rate_q18,
                          -short_short_support, -short_pre_support,
                          -short_post_support,
                          -not_long_repeat,
                          -short_short_support_sum, -short_pre_support_sum,
                          -short_post_support_sum,
                          -short_support_length_adj_sum,
                          -pre_support_length_adj_sum,
                          -post_support_length_adj_sum,
                          -half_length_adj_sum, -total_length_pre_adj_sum,
                          -total_length_post_adj_sum
  )
  return(msec)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
