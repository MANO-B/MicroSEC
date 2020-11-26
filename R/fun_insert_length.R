#' Insert length checker for BAM file
#'
#' This function attempts to check the insert length in BAM file.
#'
#' @param bam_file Path of the BAM file.
#' @return df_bam
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Biostrings DNAStringSet
#' @examples
#' fun_insert_length(
#'   system.file("extdata", "test.bam", package = "MicroSEC")
#' )
#' @export
fun_insert_length <- function(bam_file) {
  file_bam <- BamFile(bam_file)
  df_bam <- scanBam(file_bam)[[1]]
  return(stats::median(df_bam$isize[df_bam$isize > 0 &
                             df_bam$isize < 1000 &
                             !is.na(df_bam$isize)]))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
