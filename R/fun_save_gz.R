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
    data.table::fwrite(MSEC, file=OUTPUT, sep = "\t")
  } else{
    data.table::fwrite("", file=OUTPUT, sep = "\t")
  }
}


# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
