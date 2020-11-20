#' Read check function.
#'
#' This function attempts to check the mutation profile in each read.
#'
#' @param df_mutation Mutation information.
#' @param df_BAM Data from the BAM file.
#' @param df_mut_call Read ID list.
#' @param genome Reference genome for the data.
#' @param Chr_No Chromosome number.
#' @param SAMPLE_NAME Sample name (character)
#' @param READ_length The read length in the sequence.
#' @param ADAPTER_SEQ_1 The Read 1 adapter sequence of the library.
#' @param ADAPTER_SEQ_2 The Read 2 adapter sequence of the library.
#' @param Short_Homology_search_length Small sequence for homology search.
#' @param PROGRESS_BAR "Y": You can see the progress visually.
#' @return list(MSEC, Homology_search)
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom stringr str_split
#' @importFrom stringr str_sub
#' @importFrom stringr str_replace
#' @importFrom BiocGenerics which.max
#' @importFrom gtools asc
#' @importFrom Biostrings matchPattern
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings trimLRPatterns
#' @importFrom Biostrings width
#' @importFrom Biostrings start
#' @importFrom BiocGenerics as.data.frame
#' @export
fun_read_check = function(df_mutation,
                          df_BAM,
                          df_mut_call,
                          genome,
                          Chr_No,
                          SAMPLE_NAME,
                          READ_length,
                          ADAPTER_SEQ_1,
                          ADAPTER_SEQ_2,
                          Short_Homology_search_length,
                          PROGRESS_BAR){
  if(length(df_mutation[,1]) > 0){
    # initialize
    MSEC = NULL
    Homology_search = NULL
    Mut_depth = NULL
    Length_Flag = 0
    CHROM = ""
    Pre_search_length = 0
    Pre_search_length_default = 0
    Post_search_length_default = 20
    Post_search_length = 20
    Minimum_Hairpin_length = 15
    Width = 150
    Hairpin_search_length = 50
    Max_mutation_search = 50
    neighbor_length = 20
    Laxness = 1
    
    # analyze each somatic mutation
    for(i in 1:length(df_mutation[,1])){
      Error = 0
      indel_status = 0
      indel_length = 0
      WITH_INDEL_1 = FALSE
      WITH_INDEL_2 = FALSE
      mut_call = logical(0)
      Co_mut_Pre = 3
      Co_mut_Post = 3
      penalty_Pre = 0
      penalty_Post = 0
      penalty_Post = 0
      near_indel_Pre = 0
      near_indel_Pre_candidate = 0
      near_indel_Post = 0
      near_indel_Post_candidate = 0
      Mut_depth_tmp = rep(0, 160)
      Caution = ""
      # extract mutation supporting reads
      if(df_mutation[i,"Chr"] != CHROM){
        CHROM = df_mutation[i,"Chr"]
        ID_No = df_BAM$rname == CHROM
        df_BAM_qname = df_BAM$qname[ID_No]
        df_BAM_seq = df_BAM$seq[ID_No]
        df_BAM_strand = df_BAM$strand[ID_No]
        df_BAM_cigar = df_BAM$cigar[ID_No]
        df_BAM_qual = df_BAM$qual[ID_No]
        df_BAM_pos = df_BAM$pos[ID_No]
      }
      mut_read = df_mut_call %>% 
        filter(Chr == df_mutation[i,"Chr"] & Pos == df_mutation[i,"Pos"])
      mut_type = str_split(df_mutation[i,"Mut_type"],"-")[[1]][[2]]
      if(mut_type == "ins"){
        indel_length = nchar(df_mutation[i,"Alt"]) - 1
        indel_status = 1
        WITH_INDEL_1 = TRUE
        WITH_INDEL_2 = TRUE
      } else if(mut_type == "del"){
        indel_length = nchar(df_mutation[i,"Ref"]) - 1
        indel_status = 1
        WITH_INDEL_1 = TRUE
        WITH_INDEL_2 = TRUE
      }
      if(dim(mut_read)[[1]] == 0 & indel_status == 1){
        for(tmp in 1:Max_mutation_search){
          if(dim(mut_read)[[1]] == 0){
            mut_read = df_mut_call %>%
              filter(Chr == df_mutation[i,"Chr"] & 
                     Pos == (df_mutation[i,"Pos"] - tmp))
            Error = tmp
          }
        }
      }
      if(dim(mut_read)[1] > 0){
        mut_detail = str_split(mut_read$Mut, pattern = ";")[[1]]
        mut_read_ID_list = str_split(mut_read$Mut_ID, pattern = ";")[[1]]
        if(indel_status == 0){
          mut_call = which(mut_detail == str_sub(df_mutation[i,"Alt"], 
                                                 start = 1, end = 1))
        } else if(mut_type == "del"){
          mut_call = which(mut_detail ==
                             str_replace(df_mutation[i,"Ref"],
                                         pattern=df_mutation[i,"Alt"],
                                         replacement=".-"))
        } else if(mut_type == "ins"){
          mut_call = which(mut_detail ==
                             str_replace(df_mutation[i,"Alt"],
                                         pattern=df_mutation[i,"Ref"],
                                         replacement=".+"))
        }
      }
      # if mutation supporting reads exist
      if(length(mut_call) > 0){
        Pre_support_length = 0
        Post_support_length = 0
        short_support_length = 0
        Pre_Minimum_length = READ_length
        Post_Minimum_length = READ_length
        Pre_rep_status = 0
        Post_rep_status = 0
        Homopolymer_status = 0
        FLAG_Hairpin = 0
        Hairpin_length = 0
        Soft_Clipped_read = 0
        Total_read = 0
        Low_quality_base = 0
        Homology_search_tmp = NULL
        distant_homology = 0
        mutated_ID = mut_read_ID_list[mut_call]
        mut_read_ID = str_sub(str_split(mutated_ID,",")[[1]], start=1, end=-2)
        mut_read_strand = str_sub(str_split(mutated_ID,",")[[1]],
                                  start=-1, end=-1)
        mut_read_strand = gsub("r", "-", gsub("f", "+", mut_read_strand))
        Ref_seq = genome[[df_mutation[i,"Chr"]]][
          (df_mutation[i,"Pos"] - Error - Width):
          (df_mutation[i,"Pos"] - Error + Width)]
        Ref_indel = c(Ref_seq[1:Width], DNAString(df_mutation[i,"Alt"]), 
                      Ref_seq[(Width + nchar(df_mutation[i,"Ref"]) + 1):
                                (2 * Width + 1)])
        Neighbor_seq = df_mutation[i,"Neighborhood_sequence"]
        Alt_length = nchar(df_mutation[i,"Alt"])
        # sequence information around the mutation position
        Setting = fun_setting(PRE = Pre_search_length_default,
                              POST = Post_search_length_default,
                              Neighbor_seq = Neighbor_seq,
                              neighbor_length = neighbor_length,
                              Alt_length = Alt_length)
        Pre_search_length = Setting[[1]]
        Post_search_length = Setting[[2]]
        Peri_seq_1 = Setting[[3]]
        Peri_seq_2 = Setting[[4]]
        near_list_1 = df_mut_call %>%
          filter(Chr == df_mutation[i,"Chr"] &
                 Pos >= (df_mutation[i,"Pos"] - Error - Pre_search_length) &
                 Pos <= (df_mutation[i,"Pos"] + Post_search_length) &
                 Pos != (df_mutation[i,"Pos"] - Error))
        near_list_2 = df_mut_call %>%
          filter(Chr == df_mutation[i,"Chr"] &
                 Pos >= (df_mutation[i,"Pos"] - Error - Post_search_length) &
                 Pos <= (df_mutation[i,"Pos"] + Pre_search_length) &
                 Pos != (df_mutation[i,"Pos"] - Error))
        Mut_near_1 = dim(near_list_1)[1]
        Mut_near_2 = dim(near_list_2)[1]
        if(all(nchar(near_list_1$Ref) != nchar(near_list_1$Alt))){
          WITH_INDEL_1 = TRUE
        }
        if(all(nchar(near_list_2$Ref) != nchar(near_list_2$Alt))){
          WITH_INDEL_2 = TRUE
        }
        # short repeat around indel mutations
        if(mut_type == "ins"){
          Rep_status = fun_repeat_check(Rep_A = DNAString(df_mutation[i,"Ref"]),
                                        Rep_B = DNAString(df_mutation[i,"Alt"]),
                                        Ref_seq = Ref_seq,
                                        Width = Width,
                                        Del = 0)
          Pre_rep_status = Rep_status[[1]]
          Post_rep_status = Rep_status[[2]]
          Homopolymer_status = Rep_status[[3]]
        } else if(mut_type == "del"){
          Rep_status = fun_repeat_check(Rep_A = DNAString(df_mutation[i,"Alt"]),
                                        Rep_B = DNAString(df_mutation[i,"Ref"]),
                                        Ref_seq = Ref_seq,
                                        Width = Width,
                                        Del = 1)
          Pre_rep_status = Rep_status[[1]]
          Post_rep_status = Rep_status[[2]]
          Homopolymer_status = Rep_status[[3]]
        }
        for(depth in 1:(READ_length + 1)){
          Mut_depth_tmp[depth] = 
            sum(df_BAM_pos == (df_mutation[i,"Pos"] + 1 - depth))
        }
        if(PROGRESS_BAR == "Y"){
          pb = txtProgressBar(min = 0, 
                              max = max(1, length(mut_read_ID)),
                              width = 20,
                              style = 3)
          pb_t = ceiling(length(mut_read_ID) / 100)
        }
        # analyze each mutation supporting read
        for(j in 1:length(mut_read_ID)){
          if(Length_Flag == 1){
            Setting = fun_setting(PRE = Pre_search_length_default,
                                  POST = Post_search_length_default,
                                  Neighbor_seq = Neighbor_seq,
                                  neighbor_length = neighbor_length,
                                  Alt_length = Alt_length)
            Pre_search_length = Setting[[1]]
            Post_search_length = Setting[[2]]
            Peri_seq_1 = Setting[[3]]
            Peri_seq_2 = Setting[[4]]
          }
          Length_Flag = 0
          mut_position = 0
          FLAG_1 = Pre_search_length
          FLAG_2 = Pre_search_length
          # progress bar
          if(PROGRESS_BAR == "Y"){
            setTxtProgressBar(pb, j)
            if((j - 1) %% pb_t == 0){
              cat(paste(" Mutation:", i, "/",length(df_mutation[,1]), "      "))
            }
          }
          # specific read selection
          ID_No = df_BAM_qname == mut_read_ID[[j]]
          df_seq = df_BAM_seq[ID_No]
          df_strand = df_BAM_strand[ID_No]
          df_cigar = df_BAM_cigar[ID_No]
          df_qual = df_BAM_qual[ID_No]
          df_seq = df_seq[df_strand == mut_read_strand[[j]]]
          df_cigar = df_cigar[df_strand == mut_read_strand[[j]]]
          df_qual = df_qual[df_strand == mut_read_strand[[j]]]
          if(length(df_seq) > 1){
            df_cigar = df_cigar[which.max(width(df_seq))]
            df_qual = df_qual[which.max(width(df_seq)),]
            df_seq = df_seq[which.max(width(df_seq)),]
          }
          df_seq = df_seq[[1]]
          df_qual = as.vector(asc(as.character(df_qual[1])))
          if(length(df_seq) > 0){
            
            # determine mutation position in each read
            mutation_supporting_1 =
              matchPattern(Peri_seq_1,
                           df_seq,
                           max.mismatch=Mut_near_1 + Laxness,
                           min.mismatch=0, 
                           with.indels=WITH_INDEL_1,
                           fixed=FALSE)
            mutation_supporting_2 = 
              matchPattern(Peri_seq_2,
                           df_seq, 
                           max.mismatch=Mut_near_2 + Laxness,
                           min.mismatch=0,
                           with.indels=WITH_INDEL_2, 
                           fixed=FALSE)
            if(length(mutation_supporting_1) != 1 & 
               length(mutation_supporting_2) != 1){
              Length_Flag = 1
              if(length(mutation_supporting_1) == 0 & 
                 length(mutation_supporting_2) == 0){
                for(Lax_1 in seq(3, 9, length=3)){
                  for(Lax_2 in 2:0){
                    if(length(mutation_supporting_1) != 1 &
                       length(mutation_supporting_2) != 1){
                      Setting = fun_setting(
                        PRE = Pre_search_length_default,
                        POST = Post_search_length_default - Lax_1 * Laxness,
                        Neighbor_seq = Neighbor_seq,
                        neighbor_length = neighbor_length,
                        Alt_length = Alt_length)
                      Pre_search_length = Setting[[1]]
                      Post_search_length = Setting[[2]]
                      Peri_seq_1 = Setting[[3]]
                      Peri_seq_2 = Setting[[4]]
                      mutation_supporting_1 = matchPattern(
                        Peri_seq_1, df_seq,
                        max.mismatch=Mut_near_1 + Lax_2 * Laxness,
                        min.mismatch=0,
                        with.indels=WITH_INDEL_1,
                        fixed=FALSE)
                      mutation_supporting_2 = matchPattern(
                        Peri_seq_2,
                        df_seq,
                        max.mismatch=Mut_near_2 + Lax_2 * Laxness,
                        min.mismatch=0,
                        with.indels=WITH_INDEL_2, 
                        fixed=FALSE)
                    }
                  }
                }
              }
              if(length(mutation_supporting_1) != 1 &
                 length(mutation_supporting_2) != 1){
                for(Lax_1 in seq(5, 15, length=3)){
                  for(Lax_2 in 3:0){
                    if(length(mutation_supporting_1) != 1 &
                       length(mutation_supporting_2) != 1){
                      Setting = fun_setting(
                        PRE = Pre_search_length_default + Lax_1 * Laxness,
                        POST = Post_search_length_default,
                        Neighbor_seq = Neighbor_seq,
                        neighbor_length = neighbor_length,
                        Alt_length = Alt_length)
                      Pre_search_length = Setting[[1]]
                      Post_search_length = Setting[[2]]
                      Peri_seq_1 = Setting[[3]]
                      Peri_seq_2 = Setting[[4]]
                      mutation_supporting_1 = matchPattern(
                        Peri_seq_1,
                        df_seq,
                        max.mismatch=Mut_near_1 + Lax_2 * Laxness, 
                        min.mismatch=0,
                        with.indels=WITH_INDEL_1, 
                        fixed=FALSE)
                      mutation_supporting_2 = matchPattern(
                        Peri_seq_2, 
                        df_seq,
                        max.mismatch=Mut_near_2 + Lax_2 * Laxness,
                        min.mismatch=0, 
                        with.indels=WITH_INDEL_2,
                        fixed=FALSE)
                    }
                  }
                }
              }
            }
            if(length(mutation_supporting_1) == 1){
              mut_position = min(
                length(df_seq),
                start(mutation_supporting_1) + Pre_search_length)
              FLAG_1 = Hairpin_search_length
            } else if(length(mutation_supporting_2) == 1){
              mut_position = min(
                length(df_seq),
                start(mutation_supporting_2) + Post_search_length)
              FLAG_2 = Hairpin_search_length
            }
            if(mut_position > 0){
              Total_read = Total_read + 1
              Reverse_seq = logical(0)
              FLAG_Hairpin_tmp = 0
              check_hairpin = 1
              
              # search co-mutations on neighbor
              if(mut_position > 10){
                near_indel_Pre_candidate = near_indel_Pre_candidate + 1
                Co_mut_Pre_tmp = length(
                  matchPattern(
                    df_seq[(mut_position - 10):mut_position],
                    Ref_indel[(Width - 9):(Width + 1)],
                    max.mismatch=0, 
                    min.mismatch=3,
                    with.indels=FALSE, 
                    fixed=TRUE))
                if(Co_mut_Pre_tmp == 0){
                  near_indel_Pre = near_indel_Pre + 1
                }
              }
              
              if(mut_position < length(df_seq) - 10){
                near_indel_Post_candidate = near_indel_Post_candidate + 1
                Co_mut_Post_tmp = length(
                  matchPattern(
                    df_seq[mut_position:(mut_position + 10)],
                    Ref_indel[(Width + 1):(Width + 11)],
                    max.mismatch=0, 
                    min.mismatch=3,
                    with.indels=FALSE, 
                    fixed=TRUE))
                if(Co_mut_Post_tmp == 0){
                  near_indel_Post = near_indel_Post + 1
                }
              }

              if(indel_status == 1){
                comut_FLAG = TRUE
                for(comut in 0:3){
                  if(comut_FLAG == TRUE)
                    Co_mut_Pre_tmp = length(
                      matchPattern(df_seq[
                        max(1,(mut_position - indel_length - 1)):
                          min((mut_position + nchar(df_mutation[i,"Alt"]) + 4),
                              length(df_seq))],
                        Ref_indel[(Width - indel_length):
                                  (Width + nchar(df_mutation[i,"Alt"]) + 5)],
                         max.mismatch=comut, 
                         min.mismatch=comut,
                         with.indels=FALSE, 
                         fixed=TRUE))
                  if(Co_mut_Pre_tmp > 0){
                    Co_mut_Pre = min(Co_mut_Pre, comut)
                    comut_FLAG = FALSE
                  }
                }
                comut_FLAG = TRUE
                for(comut in 0:3){
                  if(comut_FLAG == TRUE)
                    Co_mut_Post_tmp = length(
                      matchPattern(df_seq[
                        max(1,(mut_position - 4)):
                        min((mut_position + nchar(df_mutation[i,"Alt"]) +
                               indel_length - 1),
                            length(df_seq))],
                        Ref_indel[(Width - 3):
                                  (Width + nchar(df_mutation[i,"Alt"]) +
                                     indel_length)],
                        max.mismatch=comut, 
                        min.mismatch=comut,
                        with.indels=FALSE, 
                        fixed=TRUE))
                  if(Co_mut_Post_tmp > 0){
                    Co_mut_Post = min(Co_mut_Post, comut)
                    comut_FLAG = FALSE
                  }
                }
              } else{
                if((mut_position - max(1, Alt_length * 4 - 5)) > 0){
                  comut_FLAG = TRUE
                  for(comut in 0:3){
                    if(comut_FLAG == TRUE)
                      Co_mut_Pre_tmp = length(
                        matchPattern(df_seq[
                          (mut_position - max(1, Alt_length * 4 - 5)):
                            mut_position],
                          Ref_indel[(Width - max(1, Alt_length * 4 - 5) + 1):
                                      (Width + 1)],
                          max.mismatch=comut, 
                          min.mismatch=comut,
                          with.indels=FALSE, 
                          fixed=TRUE))
                    if(Co_mut_Pre_tmp > 0){
                      Co_mut_Pre = min(Co_mut_Pre, comut)
                      comut_FLAG = FALSE
                    }
                  }
                  comut_FLAG = TRUE
                  for(comut in 0:3){
                    if(comut_FLAG == TRUE)
                      Co_mut_Post_tmp = length(
                        matchPattern(df_seq[
                          mut_position:
                            (mut_position + max(1, Alt_length * 4 - 5))],
                          Ref_indel[(Width + 1):
                                      (Width + max(1, Alt_length * 4 - 5) + 1)],
                          max.mismatch=comut, 
                          min.mismatch=comut,
                          with.indels=FALSE,
                          fixed=TRUE))
                    if(Co_mut_Post_tmp > 0){
                      Co_mut_Post = min(Co_mut_Post, comut)
                      comut_FLAG = FALSE
                    }
                  }
                }
              }
  
              # hairpin length calculation
              Hairpin_seq = fun_hairpin_trimming(
                          df_seq[max(1,(mut_position - FLAG_1)):
                                 min(length(df_seq),
                                     (mut_position + FLAG_2 + Alt_length - 1))],
                          mut_read_strand[[j]],
                          ADAPTER_SEQ_1, ADAPTER_SEQ_2)
              if(Minimum_Hairpin_length < length(Hairpin_seq)){
                for(hair in Minimum_Hairpin_length:length(Hairpin_seq)){
                  if(check_hairpin == 1){
                    hairpin_status = fun_hairpin_check(
                                          Hairpin_seq_tmp = Hairpin_seq[1:hair], 
                                          Ref_seq, 
                                          Hairpin_length,
                                          hair)
                    Hairpin_length = hairpin_status[[1]]
                    check_hairpin = hairpin_status[[2]]
                    FLAG_Hairpin_tmp = max(FLAG_Hairpin_tmp, check_hairpin)
                  }
                }
              }
              if(FLAG_Hairpin_tmp == 0){
                check_hairpin = 1
                Hairpin_seq = fun_hairpin_trimming(
                  df_seq[max(1,(mut_position - Minimum_Hairpin_length)):
                  min(length(df_seq),
                     (mut_position + Minimum_Hairpin_length + Alt_length - 1))],
                  mut_read_strand[[j]],
                  ADAPTER_SEQ_1, ADAPTER_SEQ_2)
                if(Minimum_Hairpin_length < length(Hairpin_seq)){
                  for(hair in Minimum_Hairpin_length:length(Hairpin_seq)){
                    if(check_hairpin == 1){
                      hairpin_status = fun_hairpin_check(
                        Hairpin_seq[1:hair], 
                        Ref_seq, 
                        Hairpin_length,
                        hair)
                      Hairpin_length = hairpin_status[[1]]
                      check_hairpin = hairpin_status[[2]]
                      hairpin_status = fun_hairpin_check(
                        Hairpin_seq[(length(Hairpin_seq) - hair + 1):
                                     length(Hairpin_seq)], 
                        Ref_seq, 
                        Hairpin_length,
                        hair)
                      Hairpin_length = hairpin_status[[1]]
                      check_hairpin = max(check_hairpin, hairpin_status[[2]])
                      FLAG_Hairpin_tmp = max(FLAG_Hairpin_tmp, check_hairpin)
                    }
                  }
                }
              }
              FLAG_Hairpin = FLAG_Hairpin + FLAG_Hairpin_tmp
              
              # supporting length calculation
              support_status = fun_support(df_cigar,
                                     df_seq,
                                     mut_read_strand[[j]],
                                     ADAPTER_SEQ_1, ADAPTER_SEQ_2,
                                     mut_position,
                                     Alt_length,
                                     indel_status)
              Pre_support_length_tmp = support_status[[1]]
              Post_support_length_tmp = support_status[[2]]
              Soft_Clipped_read_tmp = support_status[[3]]
  
              # adjustment for consecutive snv
              if(mut_type == "snv"){
                Post_support_length_tmp =
                  Post_support_length_tmp - Alt_length + 1
              }
              
              # read quality check
              Low_quality_base = 
                Low_quality_base + sum(df_qual < 51) / length(df_seq)
              
              # save sequence for homologous region search
              if(FLAG_Hairpin_tmp == 0 & FLAG_Hairpin == 0){
                if(indel_status == 1){
                  Pre_Homology_search_seq = df_seq[1:
                    min(length(df_seq), 
                        mut_position + Short_Homology_search_length +
                                       Post_rep_status + Alt_length - 1)]
                  Post_Homology_search_seq = df_seq[
                    max(1, mut_position - Short_Homology_search_length
                                        - Pre_rep_status):
                    length(df_seq)]
                } else{
                  Pre_Homology_search_seq = df_seq[1:
                    min(length(df_seq), 
                     mut_position + Short_Homology_search_length + Alt_length)]
                  Post_Homology_search_seq = df_seq[
                    max(1, mut_position - Short_Homology_search_length):
                    length(df_seq)]
                }
                if(mut_read_strand[[j]] == "-"){
                  Pre_Homology_search_seq = fun_hairpin_trimming(
                    Pre_Homology_search_seq,
                    mut_read_strand[[j]],
                    ADAPTER_SEQ_1, ADAPTER_SEQ_2)
                  Pre_Homology_search_seq = 
                    reverseComplement(Pre_Homology_search_seq)
                } else{
                  Post_Homology_search_seq = fun_hairpin_trimming(
                    Post_Homology_search_seq,
                    mut_read_strand[[j]],
                    ADAPTER_SEQ_1, ADAPTER_SEQ_2)
                  Post_Homology_search_seq =
                    reverseComplement(Post_Homology_search_seq)
                }
                Homology_search_tmp = rbind(Homology_search_tmp, 
                  data.frame(SAMPLE_NAME,
                             Chr = df_mutation[i,"Chr"],
                             Pos = df_mutation[i,"Pos"], 
                             Ref = df_mutation[i,"Ref"], 
                             Alt = df_mutation[i,"Alt"],
                             Direction = "Pre", 
                             Seq = as.character(Pre_Homology_search_seq)))
                Homology_search_tmp = rbind(Homology_search_tmp,
                  data.frame(SAMPLE_NAME, 
                             Chr = df_mutation[i,"Chr"], 
                             Pos = df_mutation[i,"Pos"], 
                             Ref = df_mutation[i,"Ref"], 
                             Alt = df_mutation[i,"Alt"], 
                             Direction = "Post", 
                             Seq = as.character(Post_Homology_search_seq)))
              }
              
              # summary
              Soft_Clipped_read = Soft_Clipped_read + Soft_Clipped_read_tmp
              Pre_support_length = max(Pre_support_length,
                                       Pre_support_length_tmp)
              Post_support_length = max(Post_support_length, 
                                        Post_support_length_tmp)
              short_support_length = 
                max(short_support_length, 
                    min(Pre_support_length_tmp, Post_support_length_tmp))
              Pre_Minimum_length = min(Pre_Minimum_length, 
                                       Pre_support_length_tmp)
              Post_Minimum_length = min(Post_Minimum_length, 
                                        Post_support_length_tmp)
            }
          }
        }
        if(indel_status == 1){
          penalty_Pre = 5 * Co_mut_Pre
          penalty_Post = 5 * Co_mut_Post
        } else{
          penalty_Pre = max(0, 4 * Alt_length - 5) + 5 * Co_mut_Pre
          penalty_Post = max(0, 4 * Alt_length - 5) + 5 * Co_mut_Post
        }

        if(fun_zero(near_indel_Pre, near_indel_Pre_candidate) == 1){
          Pre_Minimum_length = 0
          Caution = paste(Caution, 
                          "Pre_Minimum_length are set to be 0 because of neighbor co-mutations,")
          
        }
        if(fun_zero(near_indel_Post, near_indel_Post_candidate) == 1){
          Post_Minimum_length = 0
          Caution = paste(Caution, 
                          "Post_Minimum_length are set to be 0 because of neighbor co-mutations,")
        }
        
        # data formatting
        MSEC_tmp = df_mutation[i,] %>% dplyr::mutate(
          READ_length = READ_length, 
          mut_type = mut_type,
          Alt_length = Alt_length, 
          Total_read = Total_read, 
          Soft_Clipped_read = Soft_Clipped_read, 
          FLAG_Hairpin = FLAG_Hairpin, 
          Hairpin_length = Hairpin_length, 
          Pre_support_length = Pre_support_length, 
          Post_support_length = Post_support_length, 
          short_support_length = short_support_length, 
          Pre_Minimum_length = Pre_Minimum_length, 
          Post_Minimum_length = Post_Minimum_length, 
          Low_quality_base_rate_under_Q18 = fun_zero(Low_quality_base, Total_read),
          Pre_rep_status = Pre_rep_status, 
          Post_rep_status = Post_rep_status, 
          Homopolymer_status = Homopolymer_status,
          indel_status = indel_status,
          indel_length = indel_length,
          distant_homology = distant_homology,
          penalty_Pre = penalty_Pre,
          penalty_Post = penalty_Post,
          Caution = Caution
        )
        if(FLAG_Hairpin == 0){
          Homology_search = rbind(Homology_search, Homology_search_tmp)
        }
      } else{
        MSEC_tmp = df_mutation[i,] %>% dplyr::mutate(
          READ_length = READ_length, 
          mut_type = mut_type,
          Alt_length = 0, 
          Total_read = 0, 
          Soft_Clipped_read = 0, 
          FLAG_Hairpin = 0, 
          Hairpin_length = 0, 
          Pre_support_length = 0, 
          Post_support_length = 0, 
          short_support_length = 0, 
          Pre_Minimum_length = 0, 
          Post_Minimum_length = 0, 
          Low_quality_base_rate_under_Q18 = 0,
          Pre_rep_status = 0, 
          Post_rep_status = 0, 
          Homopolymer_status = 0,
          indel_status = 0,
          indel_length = 0,
          distant_homology = 0,
          penalty_Pre = 0,
          penalty_Post = 0,
          Caution = ""
        )
      }
      Mut_depth_tmp = c(0, Mut_depth_tmp)
      Mut_depth_tmp = t(cumsum(Mut_depth_tmp))
      colnames(Mut_depth_tmp) = paste("Depth", 0:160, sep="")
      MSEC = rbind(MSEC, MSEC_tmp)
      Mut_depth = rbind(Mut_depth, Mut_depth_tmp)
    }
    return(list(MSEC, Homology_search, Mut_depth))
  } else{
    return(list(df_mutation, NULL, NULL))
  }
}




# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
