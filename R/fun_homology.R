#' Homology check function.
#'
#' This function attempts to search the homologous regions.
#'
#' @param MicroSEC Mutation filtering information.
#' @param df_distant Sequences to be checked.
#' @param Minimum_Homology_search_length Minimum length to define "homologous".
#' @param PROGRESS_BAR "Y": You can see the progress visually.
#' @return MicroSEC
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom stringr str_sub
#' @importFrom stringr str_detect
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings PDict
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings countPDict
#' @importFrom BiocGenerics as.data.frame
#' @importFrom GenomeInfoDb seqnames
#' @export
fun_homology = function(MicroSEC,
                        df_distant,
                        Minimum_Homology_search_length,
                        PROGRESS_BAR){
  # initialize
  MicroSEC$distant_homology = 0
  df_distant = df_distant %>% dplyr::mutate(
    Seq = as.character(Seq))
  df_distant = df_distant %>% dplyr::mutate(
    length = nchar(Seq))
  df_distant$Number = floor(seq(1,(dim(df_distant)[1]/2 + 0.5),0.5))
  MAX_No = max(df_distant$Number)
  df_distant = df_distant[df_distant$length >= Minimum_Homology_search_length,]
  df_distant = df_distant %>% dplyr::mutate(
    fixed_seq_Pre = str_sub(
      Seq, length - Minimum_Homology_search_length + 1, length),
    fixed_seq_Post = str_sub(
      Seq, 1, Minimum_Homology_search_length)
  )
  df_distant = df_distant %>% dplyr::mutate(
    fixed_seq = ifelse(Direction == "Pre", fixed_seq_Pre, fixed_seq_Post)
  )
  df_distant = df_distant[!str_detect(df_distant$fixed_seq, pattern="N"),]
  Search_Seq_f = DNAStringSet(df_distant$fixed_seq)
  Search_Seq_r = reverseComplement(Search_Seq_f)
  Search_Seq_f = PDict(Search_Seq_f)
  Search_Seq_r = PDict(Search_Seq_r)
  distant_homology = rep(0, length(Search_Seq_f))
  if(PROGRESS_BAR == "Y"){
    pb = txtProgressBar(min = 0,
                        max = Chr_No,
                        width = 20,
                        style = 3)
  }
  for(seqname in 1:Chr_No){
    if(PROGRESS_BAR == "Y"){
      setTxtProgressBar(pb, seqname)
      cat(paste(" Chromosome screening:", seqname, "/", Chr_No, "  "))
   }
    
    target = genome[[seqnames(genome)[[seqname]]]]
    distant_homology = distant_homology + countPDict(Search_Seq_f, target)
    distant_homology = distant_homology + countPDict(Search_Seq_r, target)
  }
  df_distant$distant_homology = distant_homology
  if(PROGRESS_BAR == "Y"){
    pb = txtProgressBar(min = 0, 
                        max = max(1, dim(df_distant)[1]), 
                        width = 20,
                        style = 3)
    pb_t = ceiling(MAX_No / 100)
  }
  for(i in 1:MAX_No){
    if(PROGRESS_BAR == "Y"){
      setTxtProgressBar(pb, i)
      if((i - 1) %% pb_t == 0){
        cat(paste(" Homology count:", i, "/", MAX_No, "      "))
      }
    }
    if(dim(df_distant[df_distant$Number == i,])[1] == 1){
      tmp_distant = df_distant[df_distant$Number == i,][1,]
      if(tmp_distant$distant_homology > 0){
        MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME & 
                   MicroSEC$Chr == tmp_distant$Chr & 
                   MicroSEC$Pos == tmp_distant$Pos & 
                   MicroSEC$Ref == tmp_distant$Ref & 
                   MicroSEC$Alt == tmp_distant$Alt,]$distant_homology =  
        MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME & 
                   MicroSEC$Chr == tmp_distant$Chr & 
                   MicroSEC$Pos == tmp_distant$Pos & 
                   MicroSEC$Ref == tmp_distant$Ref & 
                   MicroSEC$Alt == tmp_distant$Alt,]$distant_homology + 1
      }
    }
    if(dim(df_distant[df_distant$Number == i,])[1] == 2){
      tmp_distant = df_distant[df_distant$Number == i,][1,]
      tmp_distant_2 = df_distant[df_distant$Number == i,][2,]
      if(tmp_distant$distant_homology > 0 | tmp_distant_2$distant_homology > 0){
        MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME & 
                   MicroSEC$Chr == tmp_distant$Chr &
                   MicroSEC$Pos == tmp_distant$Pos &
                   MicroSEC$Ref == tmp_distant$Ref & 
                   MicroSEC$Alt == tmp_distant$Alt,]$distant_homology = 
        MicroSEC[MicroSEC$Sample == tmp_distant$SAMPLE_NAME &
                   MicroSEC$Chr == tmp_distant$Chr &
                   MicroSEC$Pos == tmp_distant$Pos &
                   MicroSEC$Ref == tmp_distant$Ref & 
                   MicroSEC$Alt == tmp_distant$Alt,]$distant_homology + 1
      }
    }
  }
  return(MicroSEC)
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
