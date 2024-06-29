#' An example mutation file.
#'
#' A dataset containing the information of eight mutations.
#'
#' @format A list with 8 factors, each contains 29 variables
#' \describe{
#'   \item{Sample}{sample name}
#'   \item{Mut_type}{mutation type}
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{SimpleRepeat_TRF}{mutation locating repeat sequence}
#'   \item{Neighborhood_sequence}{neighborhood sequence}
#'   ...
#' }
"exampleMutation"

#' An example BAM file.
#'
#' A BAM file containing the information of eight mutations.
#'
#' @format A list with 8 factors, each contains 46527 variables:
#' \describe{
#'   \item{rname}{chromosome of the read}
#'   \item{qname}{read ID list}
#'   \item{seq}{sequence of the read, in DNAString}
#'   \item{strand}{strand of the read}
#'   \item{cigar}{CIGAR sequence of the read}
#'   \item{qual}{Phred quality of the read}
#'   \item{pos}{starting position of the read}
#'   \item{isize}{insert size of the read}
#'   ...
#' }
"exampleBam"

#' An example mutation file.
#'
#' A dataset containing the information of eight mutations processed by the
#' fun_read_check function.
#'
#' @format A list with 34 factors, each contains 46527 variables
#' \describe{
#'   \item{Sample}{sample name}
#'   \item{Mut_type}{mutation type}
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{SimpleRepeat_TRF}{mutation locating repeat sequence}
#'   \item{Neighborhood_sequence}{neighborhood sequence}
#'   \item{read_length}{read length}
#'   \item{mut_type}{mutation type}
#'   \item{alt_length}{length of the mutated bases}
#'   \item{total_read}{number of mutation supporting reads}
#'   \item{soft_clipped_read}{number of soft-clipped reads}
#'   \item{flag_hairpin}{number of reads produced by hairpin structure}
#'   \item{hairpin_length}{maximum length of palindromic sequences}
#'   \item{pre_support_length}{maximum 5'-supporting length}
#'   \item{post_support_length}{maximum 3'-supporting length}
#'   \item{short_support_length}{minimum supporting length}
#'   \item{pre_minimum_length}{minimum 5'-supporting length}
#'   \item{post_minimum_length}{minimum 3'-supporting length}
#'   \item{pre_minimum_length}{minimum 5'-supporting length}
#'   \item{low_quality_base_rate_under_q18}{low quality base rate}
#'   \item{low_quality_pre}{low quality base rate of 5'- side}
#'   \item{low_quality_post}{low quality base rate of 3'- side}
#'   \item{pre_farthest}{5'-farthest supported base from the mutated base}
#'   \item{post_farthest}{3'-farthest supported base from the mutated base}
#'   \item{post_rep_status}{3'-repeat sequence length}
#'   \item{homopolymer_status}{homopolymer sequence length}
#'   \item{indel_status}{whether the mutation is indel or not}
#'   \item{indel_length}{length of indel mutation}
#'   \item{distant_homology}{number of reads derived from homologous regions}
#'   \item{penalty_pre}{5'-penalty score by the mapper}
#'   \item{penalty_post}{3'-penalty score by the mapper}
#'   \item{caution}{comment}
#'   ...
#' }
"msec_read_checked"

#' An example mutation file.
#'
#' A dataset containing the information of eight mutations processed by the
#' fun_homology function.
#'
#' @format A list with 34 factors, each contains 29 variables
#' \describe{
#'   \item{Sample}{sample name}
#'   \item{Mut_type}{mutation type}
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{SimpleRepeat_TRF}{mutation locating repeat sequence}
#'   \item{Neighborhood_sequence}{neighborhood sequence}
#'   \item{read_length}{read length}
#'   \item{mut_type}{mutation type}
#'   \item{alt_length}{length of the mutated bases}
#'   \item{total_read}{number of mutation supporting reads}
#'   \item{soft_clipped_read}{number of soft-clipped reads}
#'   \item{flag_hairpin}{number of reads produced by hairpin structure}
#'   \item{hairpin_length}{maximum length of palindromic sequences}
#'   \item{pre_support_length}{maximum 5'-supporting length}
#'   \item{post_support_length}{maximum 3'-supporting length}
#'   \item{short_support_length}{minimum supporting length}
#'   \item{pre_minimum_length}{minimum 5'-supporting length}
#'   \item{post_minimum_length}{minimum 3'-supporting length}
#'   \item{pre_farthest}{5'-farthest supported base from the mutated base}
#'   \item{post_farthest}{3'-farthest supported base from the mutated base}
#'   \item{low_quality_base_rate_under_q18}{low quality base rate}
#'   \item{low_quality_pre}{low quality base rate of 5'- side}
#'   \item{low_quality_post}{low quality base rate of 3'- side}
#'   \item{pre_rep_status}{5'-repeat sequence length}
#'   \item{post_rep_status}{3'-repeat sequence length}
#'   \item{homopolymer_status}{homopolymer sequence length}
#'   \item{indel_status}{whether the mutation is indel or not}
#'   \item{indel_length}{length of indel mutation}
#'   \item{distant_homology}{number of reads derived from homologous regions}
#'   \item{penalty_pre}{5'-penalty score by the mapper}
#'   \item{penalty_post}{3'-penalty score by the mapper}
#'   \item{caution}{comment}
#'   ...
#' }
"msec_homology_searched"

#' An example mutation file.
#'
#' A dataset containing the information of eight mutations processed by the
#' fun_homology function.
#'
#' @format A list with 52 factors, each contains 29 variables
#' \describe{
#'   \item{Sample}{sample name}
#'   \item{Mut_type}{mutation type}
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{SimpleRepeat_TRF}{mutation locating repeat sequence}
#'   \item{Neighborhood_sequence}{neighborhood sequence}
#'   \item{read_length}{read length}
#'   \item{mut_type}{mutation type}
#'   \item{alt_length}{length of the mutated bases}
#'   \item{total_read}{number of mutation supporting reads}
#'   \item{soft_clipped_read}{number of soft-clipped reads}
#'   \item{flag_hairpin}{number of reads produced by hairpin structure}
#'   \item{hairpin_length}{maximum length of palindromic sequences}
#'   \item{pre_support_length}{maximum 5'-supporting length}
#'   \item{post_support_length}{maximum 3'-supporting length}
#'   \item{short_support_length}{minimum supporting length}
#'   \item{pre_minimum_length}{minimum 5'-supporting length}
#'   \item{post_minimum_length}{minimum 3'-supporting length}
#'   \item{pre_farthest}{5'-farthest supported base from the mutated base}
#'   \item{post_farthest}{3'-farthest supported base from the mutated base}
#'   \item{low_quality_base_rate_under_q18}{low quality base rate}
#'   \item{low_quality_pre}{low quality base rate of 5'- side}
#'   \item{low_quality_post}{low quality base rate of 3'- side}
#'   \item{pre_rep_status}{5'-repeat sequence length}
#'   \item{post_rep_status}{3'-repeat sequence length}
#'   \item{homopolymer_status}{homopolymer sequence length}
#'   \item{indel_status}{whether the mutation is indel or not}
#'   \item{indel_length}{length of indel mutation}
#'   \item{distant_homology}{number of reads derived from homologous regions}
#'   \item{penalty_pre}{5'-penalty score by the mapper}
#'   \item{penalty_post}{3'-penalty score by the mapper}
#'   \item{caution}{comment}
#'   \item{distant_homology_rate}{rate of reads derived from homologous regions}
#'   \item{pre_minimum_length_adj}{adjusted pre_minimum_length}
#'   \item{post_minimum_length_adj}{adjusted pre_minimum_length}
#'   \item{pre_support_length_adj}{adjusted pre_minimum_length}
#'   \item{post_support_length_adj}{adjusted pre_minimum_length}
#'   \item{shortest_support_length_adj}{the shortest short_support_length}
#'   \item{minimum_length_1}{theoretically minimum 5'-supporting length}
#'   \item{minimum_length_2}{theoretically minimum 3'-supporting length}
#'   \item{minimum_length}{theoretically minimum supporting length}
#'   \item{short_support_length_adj}{adjusted short_support_length}
#'   \item{altered_length}{substituted/inserted length}
#'   \item{half_length}{half of the read length}
#'   \item{short_support_length_total}{range of short_support_length}
#'   \item{pre_support_length_total}{range of pre_support_length}
#'   \item{post_support_length_total}{range of post_support_length}
#'   \item{half_length_total}{range of possible short_support_length}
#'   \item{total_length_total}{range of possible supporting length}
#'   \item{soft_clipped_rate}{rate of soft clipped reads}
#'   ...
#' }
"msec_summarized"

#' An example mutation file.
#'
#' A dataset containing the information of eight mutations processed by the
#' fun_homology function.
#'
#' @format A list with 37 factors, each contains 29 variables
#' \describe{
#'   \item{Sample}{sample name}
#'   \item{Mut_type}{mutation type}
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{SimpleRepeat_TRF}{mutation locating repeat sequence}
#'   \item{Neighborhood_sequence}{neighborhood sequence}
#'   \item{read_length}{read length}
#'   \item{total_read}{number of mutation supporting reads}
#'   \item{soft_clipped_read}{number of soft-clipped reads}
#'   \item{flag_hairpin}{number of reads produced by hairpin structure}
#'   \item{pre_support_length}{maximum 5'-supporting length}
#'   \item{post_support_length}{maximum 3'-supporting length}
#'   \item{short_support_length}{minimum supporting length}
#'   \item{pre_farthest}{5'-farthest supported base from the mutated base}
#'   \item{post_farthest}{3'-farthest supported base from the mutated base}
#'   \item{low_quality_base_rate_under_q18}{low quality base rate}
#'   \item{low_quality_pre}{low quality base rate of 5'- side}
#'   \item{low_quality_post}{low quality base rate of 3'- side}
#'   \item{distant_homology_rate}{rate of reads derived from homologous regions}
#'   \item{soft_clipped_rate}{rate of soft clipped reads}
#'   \item{prob_filter_1}{possibility of short-supporting length}
#'   \item{prob_filter_3_pre}{possibility of 5'-supporting length}
#'   \item{prob_filter_3_post}{possibility of 3'-supporting length}
#'   \item{filter_1_mutation_intra_hairpin_loop}{filter 1}
#'   \item{filter_2_hairpin_structure}{filter 2}
#'   \item{filter_3_microhomology_induced_mutation}{filter 3}
#'   \item{filter_4_highly_homologous_region}{filter 4}
#'   \item{filter_5_soft_clipping}{filter 5}
#'   \item{filter_6_simple_repeat}{filter 6}
#'   \item{filter_7_mutation_at_homopolymer}{filter 7}
#'   \item{filter_8_low_quality}{filter 8}
#'   \item{msec_filter_123}{any of filter 1-3}
#'   \item{msec_filter_1234}{any of filter 1-4}
#'   \item{msec_filter_all}{any of filter 1-8}
#'   \item{comment}{comment}
#'   ...
#' }
"msec_analyzed"

#' An example sequence information file.
#'
#' A dataset containing the information of reads for homology search.
#'
#' @format A list with 7 factors, each contains 1508 variables:
#' \describe{
#'   \item{sample_name}{sample name}
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{Direction}{5' (pre) or 3' (post) sequence of the mutated base}
#'   \item{Seq}{sequence for homology search}
#'   ...
#' }
"homology_searched"

#' An example sequence information file.
#'
#' A dataset containing the information of reads for homology search.
#'
#' @format Three lists with 201 factors, each contains 29 variables:
"mut_depth_checked"
