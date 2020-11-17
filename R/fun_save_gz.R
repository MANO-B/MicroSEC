#' Save function.
#'
#' This function attempts to save the filtering results.
#'
#' @param MSEC Mutation filtering information.
#' @param OUTPUT Output file name (full path).
#' @export
fun_save_gz = function(MSEC, OUTPUT){
  # save the results
  if(!is.null(MSEC)){
    gz1 <- gzfile(OUTPUT, "w")
    write.table(MSEC, gz1, sep = "\t", na="", row.names=FALSE, col.names=FALSE, quote=FALSE)
    close(gz1)
  } else{
    gz1 <- gzfile(OUTPUT, "w")
    write.table("", gz1, sep = "\t", na="", row.names=FALSE, col.names=FALSE, quote=FALSE)
    close(gz1)
  }
}


# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
