#' Insert length checker for BAM file
#'
#' This function attempts to check the insert length in BAM file.
#'
#' @param BAM_FILE Path of the BAM file.
#' @return df_BAM
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Biostrings DNAStringSet
#' @examples
#' fun_insert_length("./source/Cell_line/PC9_Cell_line_Ag_TDv4.realigned.bam")
#' @export
fun_insert_length = function(BAM_FILE){
  file_BAM = BamFile(BAM_FILE)
  df_BAM = scanBam(file_BAM)[[1]]
  return(median(df_BAM$isize[df_BAM$isize > 0 & df_BAM$isize < 1000 & !is.na(df_BAM$isize)]))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
