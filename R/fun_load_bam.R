#' BAM file loader
#'
#' This function attempts to load the BAM file.
#'
#' @param bam_file Path of the BAM file.
#' @return df_bam
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Biostrings DNAStringSet
#' @examples
#' \dontrun{
#' fun_load_bam("./source/Cell_line/PC9_Cell_line_Ag_TDv4.realigned.bam")
#' }
#' @export
fun_load_bam <- function(bam_file) {
  file_bam <- BamFile(bam_file)
  df_bam <- scanBam(file_bam)[[1]]
  df_bam <- df_bam[c("qname", "rname", "strand", "cigar", "seq", "qual", "pos")]
  return(df_bam)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
