#' Mutation data file converter
#'
#' This function attempts to convert the mutation information file.
#'
#' @param MUTATION_FILE Path of the mutation information file.
#' @param GENOME Human or Mouse genome.
#' @return df_mutation
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom BiocGenerics as.data.frame
#' @importFrom stringr str_sub
#' @importFrom stringr str_detect
#' @examples
#' fun_convert("/mnt/source/mutation.xlsx", "hg38")
#' @export
fun_convert = function(MUTATION_FILE,
                       GENOME){
  # load somatic mutation list
  df_mutation = read.xlsx(MUTATION_FILE, sheet = 1)
  # load genomic sequence
  genome = fun_load_genome(GENOME)
  # genome sequence extraction functions
  fun_genome = function(x, y){
    r = NULL
    for(i in 1:length(x)){
      r = c(r, as.character(genome[[x[i]]][y[i]]))
    }
    return(r)
  }
  fun_genome_2 = function(x, y, z){
    r = NULL
    for(i in 1:length(x)){
      r = c(r, as.character(genome[[x[i]]][y[i]:z[i]]))
    }
    return(r)
  }
  # position conversion
  df_mutation = df_mutation %>% mutate(
    Gene=Hugo_Symbol,
    HGVS.p=Protein_Change, 
    Pos=Start_Position, 
    Start=Start_Position, 
    End=End_Position,
    Mut_type=Variant_Type,
    Ref=Reference,
    Alt=Tumor_Seq)
  df_mutation = df_mutation %>% dplyr::mutate(
    Mut_type = tolower(Mut_type))
  df_mutation = df_mutation %>% dplyr::mutate(
    PRE_del = fun_genome(Chr, as.integer(Start)-1),
    PRE_ins = fun_genome(Chr, as.integer(Start)),
    POST_ins = fun_genome(Chr, as.integer(End)),
    Alt_length_1 = nchar(Ref),
    Alt_length_2 = nchar(Alt))
  df_mutation = df_mutation %>% dplyr::mutate(
    Mut_type = ifelse(Mut_type == "snp", "snv", Mut_type))
  df_mutation = df_mutation %>% dplyr::mutate(
    Alt_length = (((Alt_length_1 - Alt_length_2) + abs(Alt_length_1 - Alt_length_2)) / 2) + Alt_length_2,
    Ref_ins = ifelse(Mut_type == "ins", PRE_ins, ""),
    Ref_del = ifelse(Mut_type == "del", paste(PRE_del, Ref, sep=""), ""),
    Ref_snv = ifelse(Mut_type == "snv", Ref, ""),
    Alt_ins = ifelse(Mut_type == "ins", paste(PRE_ins, Alt, sep=""), ""),
    Alt_del = ifelse(Mut_type == "del", PRE_del, ""),
    Alt_snv = ifelse(Mut_type == "snv", Alt, ""))
  df_mutation = df_mutation %>% dplyr::mutate(
    Alt_indel = paste(Alt_ins, Alt_del, Alt_snv, sep=""),
    Ref_indel = paste(Ref_ins, Ref_del, Ref_snv, sep=""))
  df_mutation = df_mutation %>% dplyr::mutate(
    Neighbor_start_1 = as.integer(Start) - 20,
    Neighbor_end_1 = as.integer(Start) - 1,
    Neighbor_start_2 = as.integer(End) + 1,
    Neighbor_end_2 = as.integer(End) + 20)
  df_mutation = df_mutation %>% dplyr::mutate(
    Pre_Neighbor = fun_genome_2(Chr, Neighbor_start_1, Neighbor_end_1),
    Post_Neighbor = fun_genome_2(Chr, Neighbor_start_2, Neighbor_end_2))
  df_mutation = df_mutation %>% dplyr::mutate(
    Neighborhood_sequence = 
      ifelse(
        Mut_type == "ins", 
          paste(Pre_Neighbor, Alt_indel, POST_ins, 
                str_sub(Post_Neighbor, 1, 19), sep=""), 
        ifelse(
          Mut_type == "del", 
            paste(Pre_Neighbor, Post_Neighbor, sep=""),
            paste(Pre_Neighbor, Alt, Post_Neighbor, sep=""))))
  df_mutation = df_mutation %>% dplyr::mutate(
    Mut_type = paste(Alt_length, "-", Mut_type, sep=""),
    Ref = Ref_indel,
    Alt = Alt_indel,
    Pos = ifelse(str_detect(Mut_type, pattern="del"), Pos - 1, Pos)
    )
  df_mutation = df_mutation %>% dplyr::select(
    -PRE_del, -PRE_ins, -POST_ins, -Alt_length_1, -Alt_length_2,
    -Alt_length, -Ref_ins, -Ref_del, -Ref_snv, -Alt_ins, -Alt_del, -Alt_snv,
    -Alt_indel, -Ref_indel, -Neighbor_start_1, -Neighbor_start_2,
    -Neighbor_end_1, -Neighbor_end_2,
    -Pre_Neighbor, -Post_Neighbor,
    -Hugo_Symbol, -Start_Position, -End_Position, -Variant_Type,
    -Reference, -Tumor_Seq, -Protein_Change, -Start, -End)
  return(df_mutation)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL