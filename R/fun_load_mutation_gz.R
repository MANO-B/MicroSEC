#' Mutation data file loader
#'
#' This function attempts to load the mutation information file.
#'
#' @param mutation_file Path of the mutation information file.
#' @return df_mutation
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @examples
#' fun_load_mutation_gz(
#'   system.file("extdata", "test_mutation.tsv", package = "MicroSEC")
#' )
#' @export
fun_load_mutation_gz <- function(mutation_file) {
  Ref <- NULL
  Alt <- NULL
  Neighborhood_sequence <- NULL
  
  # load somatic mutation list
  df_mutation <- utils::read.csv(mutation_file,
                          stringsAsFactors = FALSE,
                          header = TRUE,
                          check.names = F,
                          sep = "\t")
  # data formatting
  if (dim(df_mutation)[1] > 1) {
    if (!"SimpleRepeat_TRF" %in% colnames(df_mutation)) {
      df_mutation$SimpleRepeat_TRF <- "NA"
    }
    if (!"Transition" %in% colnames(df_mutation)) {
      df_mutation$Transition <- "NA"
    }
    df_mutation$Pos <- as.integer(df_mutation$Pos)
    df_mutation <- df_mutation %>%
      mutate(Ref = toupper(Ref)) %>%
      mutate(Alt = toupper(Alt)) %>%
      mutate(Neighborhood_sequence = toupper(Neighborhood_sequence))
  } else {
    if (!"SimpleRepeat_TRF" %in% colnames(df_mutation)) {
      df_mutation$SimpleRepeat_TRF <- character(0)
    }
    if (!"Transition" %in% colnames(df_mutation)) {
      df_mutation$Transition <- character(0)
    }
  }
  return(df_mutation)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
