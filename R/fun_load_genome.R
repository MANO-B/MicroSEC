#' Genome loading function.
#'
#' This function attempts to load the appropriate genome.
#'
#' @param GENOME Human or Mouse genome.
#' @return genome
#' @examples
#' fun_load_genome("Human")
#' @export
fun_load_genome = function(GENOME){
  if(GENOME %in% c("Human", "hg38")){
    return(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  } else if(GENOME %in% c("Mouse", "mm10")){
    return(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  } else if(GENOME %in% c("hg19")){
    return(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  }
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
