#' Mutation data file loader
#'
#' This function attempts to load the mutation information file.
#'
#' @param MUTATION_FILE Path of the mutation information file.
#' @param SAMPLE_NAME Sample name.
#' @return df_mutation
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @examples
#' fun_load_mutation_gz("./source/CCLE_SM.gz")
#' @export
fun_load_mutation_gz = function(MUTATION_FILE){
  # load somatic mutation list
  df_mutation = data.table::fread(MUTATION_FILE, stringsAsFactors=FALSE, header=TRUE, sep="\t")
  # data formatting
  if(!"SimpleRepeat_TRF" %in% colnames(df_mutation)){
    df_mutation$SimpleRepeat_TRF = "NA"
  }
  if(!"Transition" %in% colnames(df_mutation)){
    df_mutation$Transition = "NA"
  }
  df_mutation$Pos = as.integer(df_mutation$Pos)
  df_mutation = df_mutation %>%
    filter(Sample == SAMPLE_NAME) %>%
    mutate(Ref = toupper(Ref)) %>%
    mutate(Alt = toupper(Alt)) %>% 
    mutate(Neighborhood_sequence = toupper(Neighborhood_sequence))
  return(df_mutation)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
