#' Save function.
#'
#' This function attempts to save the filtering results.
#'
#' @param MSEC Mutation filtering information.
#' @param OUTPUT Output file name (full path).
#' @export
fun_save_gz = function(MSEC, OUTPUT){
  # save the results
  gz1 <- gzfile(OUTPUT, "w")
  write.table(MSEC, gz1, sep = "\t", na="", row.names=FALSE, col.names=TRUE, quote=FALSE)
  close(gz1)
}


# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
