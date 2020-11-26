#' Mutation-supporting read ID files loader
#'
#' This function attempts to load the read ID information files.
#'
#' @param read_list Path of the read ID information
#'   directory.
#' @return df_mut_call
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @examples
#' fun_load_id(
#'   system.file("extdata", package = "MicroSEC")
#' )
#' @export
fun_load_id <- function(read_list) {
  Chr <- NULL
  Pos <- NULL
  Ref <- NULL
  Alt <- NULL
  Mut_ID <- NULL
  Mut <- NULL

  # load mutation supporting read list
  files <- data.frame(file = list.files(read_list,
                                        pattern = ".gz",
                                        full.names = T))
  df_mut_call <- NULL
  for (i in files$file) {
    df_mut_call <- rbind(df_mut_call,
                         utils::read.csv(as.character(i),
                                  header = TRUE,
                                  stringsAsFactors = FALSE,
                                  sep = "\t"))
  }
  # data formatting
  df_mut_call <- df_mut_call %>%
    select(Chr, Pos, Ref, Alt, Mut_ID, Mut) %>%
    mutate(Ref = toupper(Ref)) %>%
    mutate(Alt = toupper(Alt)) %>%
    mutate(Mut = toupper(Mut))
  df_mut_call$Pos <- as.integer(df_mut_call$Pos)
  return(df_mut_call)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
