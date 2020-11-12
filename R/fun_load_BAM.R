#' BAM file loader
#'
#' This function attempts to load the BAM file.
#'
#' @param BAM_FILE Path of the BAM file.
#' @return df_BAM
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Biostrings DNAStringSet
#' @examples
#' fun_load_BAM("./source/Cell_line/PC9_Cell_line_Ag_TDv4.realigned.bam")
#' @export
fun_load_BAM = function(BAM_FILE){
  p = ScanBamParam(what=c("qname", "rname", "strand", "cigar", "seq", "qual", "pos"))
  file_BAM = BamFile(BAM_FILE)
  df_BAM = scanBam(file_BAM, param=p)[[1]]
  df_BAM = df_BAM[c("qname", "rname", "strand", "cigar", "seq", "qual", "pos")]
  df_BAM$rname = as.character(df_BAM$rname)
  df_BAM$qname = as.character(df_BAM$qname)
  df_BAM$seq = DNAStringSet(as.character(df_BAM$seq))
  df_BAM$strand = as.character(df_BAM$strand)
  df_BAM$cigar = as.character(df_BAM$cigar)
  df_BAM$qual = as.character(df_BAM$qual)
  df_BAM$pos = as.integer(df_BAM$pos)
  return(df_BAM)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
