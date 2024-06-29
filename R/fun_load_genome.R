#' Genome loading function.
#'
#' This function attempts to load the appropriate genome.
#'
#' @param organism Human or Mouse genome.
#' @return ref_genome
#' @examples
#' fun_load_genome("Human")
#' @export
fun_load_genome <- function(organism) {
  if (organism %in% c("Human", "hg38")) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("Install BSgenome.Hsapiens.UCSC.hg38 package.",
           call. = FALSE)
    }
    ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    } else if (organism %in% c("Mouse", "mm10")) {
    if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
      stop("Install BSgenome.Mmusculus.UCSC.mm10 package.",
           call. = FALSE)
    }
    ref_genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  } else if (organism %in% c("hg19")) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
      stop("Install BSgenome.Hsapiens.UCSC.hg19 package.",
           call. = FALSE)
    }
    ref_genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }
  return(ref_genome)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
