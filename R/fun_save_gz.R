#' Save function.
#'
#' This function attempts to save the filtering results.
#'
#' @param msec Mutation filtering information.
#' @param output output file name (full path).
#' @examples
#' \donttest{
#' fun_save_gz(msec_analyzed, "./test_data.tsv")
#' }
#' @export
fun_save_gz <- function(msec, output) {
  # save the results
  if(length(grep("/", output)) > 0){
    wd = gsub(
           paste("/",
                 strsplit(output, "/")[[1]][length(strsplit(output, "/")[[1]])],
                 sep = ""),
           "",
           output)
    dir.create(file.path(wd), showWarnings = FALSE)
  }
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
