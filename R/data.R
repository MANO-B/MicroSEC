#' An example mutation file.
#'
#' A dataset containing the information of four mutations.
#'
#' @format A list with 14 variables:
#' \describe{
#'   \item{Sample}{sample name}
#'   \item{Gene}{altered gene}
#'   \item{HGVS.c}{base change}
#'   \item{HGVS.p}{protein change}
#'   \item{Mut_type}{mutation type}
#'   \item{Total_QV>=20}{total reads}
#'   \item{%Alt}{altered base ratio}
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{SimpleRepeat_TRF}{mutation locating repeat sequence}
#'   \item{Neighborhood_sequence}{neighborhood sequence}
#'   \item{Transition}{base change type}
#'   ...
#' }
"exampleMutation"

#' An example BAM file.
#'
#' A dataset containing the information of four mutations.
#'
#' @format A list with 7 factors, each contains 184 variables:
#' \describe{
#'   \item{rname}{chromosome of the read}
#'   \item{qname}{read ID list}
#'   \item{seq}{sequence of the read, in DNAString}
#'   \item{strand}{strand of the read}
#'   \item{cigar}{CIGAR sequence of the read}
#'   \item{qual}{Phred quality of the read}
#'   \item{pos}{starting position of the read}
#'   ...
#' }
"exampleBAM"

#' An example mutated read ID file.
#'
#' A dataset containing the information of mutated read.
#'
#' @format A list with 7 factors, each contains 184 variables:
#' \describe{
#'   \item{Chr}{altered chromosome}
#'   \item{Pos}{altered position}
#'   \item{Ref}{reference base}
#'   \item{Alt}{altered base}
#'   \item{Mut_ID}{mutated read ID}
#'   \item{Mut}{mutation detail}
#'   ...
#' }
"exampleMutCall"
