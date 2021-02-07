#' Mutation data file loader
#'
#' This function attempts to load the mutation information file.
#'
#' @param mutation_file Path of the mutation information file.
#' @param sample_name Sample name.
#' @return df_mutation
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom BiocGenerics as.data.frame
#' @importFrom stringr str_sub
#' @importFrom stringr str_detect
#' @examples
#' fun_load_mutation(
#'   system.file("extdata", "test_mutation.xlsx", package = "MicroSEC"),
#'  "H15-11943-1-T_TDv3"
#' )
#' @export
fun_load_mutation <- function(mutation_file,
                             sample_name) {
  Sample <- NULL
  Gene <- NULL
  `HGVS.c` <- NULL
  `HGVS.p` <- NULL
  Mut_type <- NULL
  `Total_QV>=20` <- NULL
  `%Alt` <- NULL
  Chr <- NULL
  Pos <- NULL
  Ref <- NULL
  Alt <- NULL
  SimpleRepeat_TRF <- NULL
  Neighborhood_sequence <- NULL
  Transition <- NULL
  Neighbor_start_1 <- NULL
  Neighbor_start_2 <- NULL
  eighbor_end_1 <- NULL
  eighbor_end_2 <- NULL
  Pre_Neighbor <- NULL
  Post_Neighbor <- NULL
  Mut_len <- NULL


  # load somatic mutation list
  df_mutation <<- read.xlsx(mutation_file, sheet = 1)
  # data formatting
  if (!"HGVS.c" %in% colnames(df_mutation)) {
    df_mutation$HGVS.c <<- "NA"
  }
  if (!"HGVS.p" %in% colnames(df_mutation)) {
    df_mutation$HGVS.p <<- "NA"
  }
  if (!"HGVS.c" %in% colnames(df_mutation)) {
    df_mutation$`Total_QV>=20` <<- "NA"
  }
  if (!"Total_QV>=20" %in% colnames(df_mutation)) {
    df_mutation$`Total_QV>=20` <<- "NA"
  }
  if (!"%Alt" %in% colnames(df_mutation)) {
    df_mutation$`%Alt` <<- "NA"
  }
  if (!"SimpleRepeat_TRF" %in% colnames(df_mutation)) {
    df_mutation$SimpleRepeat_TRF <<- "NA"
  }
  if (!"Transition" %in% colnames(df_mutation)) {
    df_mutation$Transition <<- "NA"
  }
  if (df_mutation$Neighborhood_sequence[[1]] == "-") {
    fun_genome <- function(x, y, z) {
      r <- NULL
      for (i in seq_len(length(x))) {
        r <- c(r, as.character(ref_genome[[x[i]]][y[i]:z[i]]))
      }
      return(r)
    }
    df_mutation$Chr <<- as.character(df_mutation$Chr)
    df_mutation$Pos <<- as.integer(df_mutation$Pos)
    mut_len <- as.integer(str_split(df_mutation[,"Mut_type"], "-",
                                    simplify = TRUE)[,1])
    mut_type <- str_split(df_mutation[,"Mut_type"], "-",
                                    simplify = TRUE)[,2]
    df_mutation <<- df_mutation %>% dplyr::mutate(
      Mut_len = mut_len,
      Mut_type = mut_type)
    df_mutation <<- df_mutation %>% dplyr::mutate(
      Neighbor_start_1 =
        ifelse(
          Mut_type == "ins", Pos - 19,
          ifelse(
            Mut_type == "del", Pos - 19, Pos - 20)),
      Neighbor_end_1 = Pos - 1,
      Neighbor_start_2 =
        ifelse(
          Mut_type == "ins", Pos + 1,
          ifelse(
            Mut_type == "del", Pos + Mut_len + 1, Pos + Mut_len)),
      Neighbor_end_2 =
        ifelse(
          Mut_type == "ins", Pos + 20,
          ifelse(
            Mut_type == "del", Pos + Mut_len + 20, Pos + Mut_len + 19)))
      df_mutation <<- df_mutation %>% dplyr::mutate(
      Pre_Neighbor = fun_genome(Chr, Neighbor_start_1, Neighbor_end_1),
      Post_Neighbor = fun_genome_2(Chr, Neighbor_start_2, Neighbor_end_2))
    df_mutation <<- df_mutation %>% dplyr::mutate(
      Neighborhood_sequence =
        paste(Pre_Neighbor, Alt, Post_Neighbor, sep = ""))
  }
  #df_mutation <<- df_mutation %>%
  #  select(Sample, Gene, HGVS.c, HGVS.p, Mut_type, `Total_QV>=20`, `%Alt`, Chr,
  #         Pos, Ref, Alt, SimpleRepeat_TRF, Neighborhood_sequence, Transition)
  df_mutation$Pos <<- as.integer(df_mutation$Pos)
  df_mutation <<- df_mutation %>%
    filter(Sample == sample_name) %>%
    mutate(Ref = toupper(Ref)) %>%
    mutate(Alt = toupper(Alt)) %>%
    mutate(Neighborhood_sequence = toupper(Neighborhood_sequence)) %>%
    arrange(Chr)
  
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
