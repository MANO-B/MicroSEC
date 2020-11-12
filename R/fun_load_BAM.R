#' BAM file loader
#'
#' This function attempts to load the BAM file.
#'
#' @param BAM_FILE Path of the BAM file.
#' @return df_BAM
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Biostrings DNAStringSet
#' @examples
#' fun_load_BAM("./source/Cell_line/PC9_Cell_line_Ag_TDv4.realigned.bam")
#' @export
fun_load_BAM = function(BAM_FILE){
  file_BAM = BamFile(BAM_FILE)
  df_BAM = scanBam(file_BAM)[[1]]
  df_BAM = df_BAM[c("qname", "rname", "strand", "cigar", "seq", "qual", "pos")]
  return(df_BAM)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
