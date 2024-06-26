% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{msec_homology_searched}
\alias{msec_homology_searched}
\title{An example mutation file.}
\format{
A list with 34 factors, each contains 29 variables
\describe{
\item{Sample}{sample name}
\item{Mut_type}{mutation type}
\item{Chr}{altered chromosome}
\item{Pos}{altered position}
\item{Ref}{reference base}
\item{Alt}{altered base}
\item{SimpleRepeat_TRF}{mutation locating repeat sequence}
\item{Neighborhood_sequence}{neighborhood sequence}
\item{read_length}{read length}
\item{mut_type}{mutation type}
\item{alt_length}{length of the mutated bases}
\item{total_read}{number of mutation supporting reads}
\item{soft_clipped_read}{number of soft-clipped reads}
\item{flag_hairpin}{number of reads produced by hairpin structure}
\item{hairpin_length}{maximum length of palindromic sequences}
\item{pre_support_length}{maximum 5'-supporting length}
\item{post_support_length}{maximum 3'-supporting length}
\item{short_support_length}{minimum supporting length}
\item{pre_minimum_length}{minimum 5'-supporting length}
\item{post_minimum_length}{minimum 3'-supporting length}
\item{pre_farthest}{5'-farthest supported base from the mutated base}
\item{post_farthest}{3'-farthest supported base from the mutated base}
\item{low_quality_base_rate_under_q18}{low quality base rate}
\item{low_quality_pre}{low quality base rate of 5'- side}
\item{low_quality_post}{low quality base rate of 3'- side}
\item{pre_rep_status}{5'-repeat sequence length}
\item{post_rep_status}{3'-repeat sequence length}
\item{homopolymer_status}{homopolymer sequence length}
\item{indel_status}{whether the mutation is indel or not}
\item{indel_length}{length of indel mutation}
\item{distant_homology}{number of reads derived from homologous regions}
\item{penalty_pre}{5'-penalty score by the mapper}
\item{penalty_post}{3'-penalty score by the mapper}
\item{caution}{comment}
...
}
}
\usage{
msec_homology_searched
}
\description{
A dataset containing the information of eight mutations processed by the
fun_homology function.
}
\keyword{datasets}
