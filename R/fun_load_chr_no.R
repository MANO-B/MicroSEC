#' Chromosome number loading function.
#'
#' This function attempts to load the chromosome number.
#'
#' @param GENOME Human or Mouse genome.
#' @return Chr_No
#' @examples
#' fun_load_chr_no("Human")
#' @export
fun_load_chr_no = function(GENOME){
  if(GENOME %in% c("Human", "hg19", "hg38")){
    return(24)
  }
  if(GENOME %in% c("Mouse", "mm10")){
    return(21)
  }
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
