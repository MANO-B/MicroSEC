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
#' fun_load_mutation("./source/CCLE.xlsx")
#' @export
fun_load_mutation = function(MUTATION_FILE,
                             SAMPLE_NAME){
  # load somatic mutation list
  df_mutation = read.xlsx(MUTATION_FILE, sheet = 1)
  # data formatting
  if(!"HGVS.c" %in% colnames(df_mutation)){
    df_mutation$HGVS.c = "NA"
  }
  if(!"HGVS.p" %in% colnames(df_mutation)){
    df_mutation$HGVS.p = "NA"
  }
  if(!"HGVS.c" %in% colnames(df_mutation)){
    df_mutation$`Total_QV>=20` = "NA"
  }
  if(!"%Alt" %in% colnames(df_mutation)){
    df_mutation$`%Alt` = "NA"
  }
  if(!"SimpleRepeat_TRF" %in% colnames(df_mutation)){
    df_mutation$SimpleRepeat_TRF = "NA"
  }
  if(!"Transition" %in% colnames(df_mutation)){
    df_mutation$Transition = "NA"
  }
  df_mutation = df_mutation %>%
    select(Sample, Gene, HGVS.c, HGVS.p, Mut_type, `Total_QV>=20`, `%Alt`, Chr, Pos, Ref, Alt, SimpleRepeat_TRF, Neighborhood_sequence, Transition)
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
