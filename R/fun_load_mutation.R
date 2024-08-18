#' Mutation data file loader
#'
#' This function attempts to load the mutation information file.
#'
#' @param mutation_file Path of the mutation information file.
#' @param sample_name Sample name.
#' @param ref_genome Reference genome for the data.
#' @param chr_no Reference genome chromosome number (human=24, mouse=22).
#' @param simple_repeat_list Optional, set simple repeat bed file path.
#' @return df_mutation
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importFrom dplyr arrange
#' @importFrom dplyr case_when
#' @importFrom BiocGenerics as.data.frame
#' @importFrom stringr str_sub
#' @importFrom stringr str_detect
#' @importFrom stats complete.cases
#' @importFrom utils read.csv
#' @examples
#' fun_load_mutation(
#'   system.file("extdata", "mutation_list.tsv", package = "MicroSEC"),
#'   "sample",
#'   BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#'   24
#' )
#' @export
fun_load_mutation <- function(mutation_file,
                             sample_name,
                             ref_genome,
                             chr_no,
                             simple_repeat_list = "") {
  Sample <- NULL
  Mut_type <- NULL
  Chr <- NULL
  Pos <- NULL
  Ref <- NULL
  Alt <- NULL
  RefLen <- NULL
  AltLen <- NULL
  Mut_len <- NULL
  SimpleRepeat_TRF <- NULL
  Neighborhood_sequence <- NULL
  Neighbor_start_1 <- NULL
  Neighbor_start_2 <- NULL
  Neighbor_end_1 <- NULL
  Neighbor_end_2 <- NULL
  Pre_Neighbor <- NULL
  Post_Neighbor <- NULL
  V1 <- NULL
  V2 <- NULL
  V3 <- NULL
  
  if (ref_genome@user_seqnames[[1]] == "chr1") {
    chromosomes <- paste0("chr", c(seq_len(chr_no - 2),"X", "Y"))
  }
  if (ref_genome@user_seqnames[[1]] == "1") {
    chromosomes <- paste0("", c(seq_len(chr_no - 2),"X", "Y"))
  }
  
  # load somatic mutation list
  df_mutation <- read.csv(mutation_file,
                          stringsAsFactors = FALSE,
                          header = TRUE,
                          check.names = F,
                          sep = "\t")
  df_mutation <- df_mutation[complete.cases(df_mutation$Sample),] %>%
    filter(Sample == sample_name)
  df_mutation <- df_mutation[order(df_mutation$Chr, df_mutation$Pos),]
  
  # data formatting
  df_mutation <- df_mutation %>% mutate(
    Alt = dplyr::case_when(
      Alt == "*" ~ str_sub(Ref, 1, 1),
      TRUE ~ Alt
    )
  )
  if (length(df_mutation$Neighborhood_sequence) > 0) {
    if (df_mutation$Neighborhood_sequence[[1]] == "-") {
      df_mutation$RefLen = nchar(df_mutation$Ref)
      df_mutation$AltLen = nchar(df_mutation$Alt)
      df_mutation = df_mutation %>% mutate(
        Mut_type = case_when(
          RefLen == AltLen ~ paste0(RefLen, "-snv"),
          RefLen > AltLen ~ paste0(RefLen-1, "-del"),
          RefLen < AltLen ~ paste0(AltLen-1, "-ins"),
        )
      ) %>% select(-RefLen, -AltLen)
      
    }
    if (!"SimpleRepeat_TRF" %in% colnames(df_mutation)) {
      df_mutation$SimpleRepeat_TRF <- "-"
    }
    
    if (file.exists(simple_repeat_list) &
        !df_mutation$SimpleRepeat_TRF[1] %in% c("Y", "N")) {
      simple_repeat_info <- read.csv(simple_repeat_list,
                                     stringsAsFactors = FALSE,
                                     header = FALSE,
                                     check.names = FALSE,
                                     sep = "\t")[,1:3]
      simple_repeat_info <- simple_repeat_info %>%
        filter(V1 %in% chromosomes)
      chr_now <- chromosomes[1]
      simple_repeat_now <- simple_repeat_info %>%
        filter(V1 == chr_now)
      for (k in seq_len(length(df_mutation$SimpleRepeat_TRF))) {
        if (df_mutation$Chr[k] != chr_now) {
          chr_now <- df_mutation$Chr[k]
          simple_repeat_now <- simple_repeat_info %>%
            filter(V1 == chr_now)
        }
        if (length((simple_repeat_now %>%
                    filter(V1 == df_mutation$Chr[k] &
                           V2 <= df_mutation$Pos[k] &
                           V3 >= df_mutation$Pos[k]))$V1)) {
          df_mutation$SimpleRepeat_TRF[k] <- "Y"
        }
      }
    }
    df_mutation = df_mutation %>% mutate(
      SimpleRepeat_TRF = case_when(
        is.na(SimpleRepeat_TRF) ~ "N",
        SimpleRepeat_TRF == "-" ~ "N",
        SimpleRepeat_TRF == "" ~ "N",
        SimpleRepeat_TRF == "Y" ~ "Y",
        SimpleRepeat_TRF == "N" ~ "N",
        SimpleRepeat_TRF == TRUE ~ "Y",
        TRUE ~ "N"
      )
    )
  
    if (df_mutation$Neighborhood_sequence[[1]] == "-") {
      fun_genome <- function(x, y, z) {
        r <- NULL
        for (i in seq_len(length(x))) {
          r <- c(r, as.character(ref_genome[[x[i]]][y[i]:z[i]]))
        }
        return(r)
      }
      df_mutation$Chr <- as.character(df_mutation$Chr)
      df_mutation$Chr_original <- df_mutation$Chr
      if (str_sub(df_mutation$Chr[[1]], start=1, end=3) != "chr" &
           ref_genome@user_seqnames[[1]] == "chr1") {
        df_mutation$Chr <- paste("chr", df_mutation$Chr, sep="")
      }
      df_mutation$Pos <- as.integer(df_mutation$Pos)
      mut_len <- as.integer(str_split(df_mutation[,"Mut_type"], "-",
                                      simplify = TRUE)[,1])
      mut_type <- str_split(df_mutation[,"Mut_type"], "-",
                                      simplify = TRUE)[,2]
      df_mutation <- df_mutation %>% mutate(
        mut_len = mut_len,
        mut_type = mut_type)
      df_mutation <- df_mutation %>% mutate(
        Neighbor_start_1 =
          ifelse(
            mut_type == "ins", Pos - 19,
            ifelse(
              mut_type == "del", Pos - 19, Pos - 20)),
        Neighbor_end_1 = Pos - 1,
        Neighbor_start_2 =
          ifelse(
            mut_type == "ins", Pos + 1,
            ifelse(
              mut_type == "del", Pos + mut_len + 1, Pos + mut_len)),
        Neighbor_end_2 =
          ifelse(
            mut_type == "ins", Pos + 20,
            ifelse(
              mut_type == "del", Pos + mut_len + 20, Pos + mut_len + 19)))
        df_mutation <- df_mutation %>% mutate(
        Pre_Neighbor = fun_genome(Chr, Neighbor_start_1, Neighbor_end_1),
        Post_Neighbor = fun_genome(Chr, Neighbor_start_2, Neighbor_end_2))
      df_mutation <- df_mutation %>% mutate(
        Neighborhood_sequence =
          paste0(Pre_Neighbor, Alt, Post_Neighbor)) %>%
        select(-Pre_Neighbor, -Post_Neighbor, -mut_len, -mut_type,
               -Neighbor_start_1, -Neighbor_start_2,
               -Neighbor_end_1, -Neighbor_end_2)
    }
    df_mutation$Pos <- as.integer(df_mutation$Pos)
    df_mutation <- df_mutation %>%
      mutate(Ref = toupper(Ref)) %>%
      mutate(Alt = toupper(Alt)) %>%
      mutate(Neighborhood_sequence = toupper(Neighborhood_sequence)) %>%
      arrange(Chr)
    return(df_mutation)
  } else {
    print("No mutation in the mutation file, exit.")
    if(NA)print("ERROR")
  }
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
