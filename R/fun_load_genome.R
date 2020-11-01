#' Genome loading function.
#'
#' This function attempts to load the appropriate genome.
#'
#' @param GENOME Human or Mouse genome.
#' @return genome
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import BSgenome.Mmusculus.UCSC.mm10
#' @examples
#' fun_load_genome("Human")
#' @export
fun_load_genome = function(GENOME){
  if(GENOME == "Human"){
    return(BSgenome.Hsapiens.UCSC.hg38)
  }
  if(GENOME == "Mouse"){
    return(BSgenome.Mmusculus.UCSC.mm10)
  }
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
