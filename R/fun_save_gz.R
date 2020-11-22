#' Save function.
#'
#' This function attempts to save the filtering results.
#'
#' @param msec Mutation filtering information.
#' @param output output file name (full path).
#' @examples
#' \dontrun{
#' fun_save_gz(msec, output)
#' }
#' @export
fun_save_gz <- function(msec, output) {
  # save the results
  gz1 <- gzfile(output, "w")
  utils::write.table(msec, gz1, sep = "\t", na = "",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  close(gz1)
}


# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
