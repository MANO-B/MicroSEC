#' Read check function.
#'
#' This function attempts to check the mutation profile in each read.
#'
#' @param df_mutation Mutation information.
#' @param df_bam Data from the BAM file.
#' @param df_mut_call Read ID list.
#' @param ref_genome Reference genome for the data.
#' @param sample_name Sample name (character)
#' @param read_length The read length in the sequence.
#' @param adapter_1 The Read 1 adapter sequence of the library.
#' @param adapter_2 The Read 2 adapter sequence of the library.
#' @param short_homology_search_length Small sequence for homology search.
#' @param progress_bar "Y": You can see the progress visually.
#' @return list(msec, homology_search)
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
#' @examples
#' \donttest{
#' fun_read_check(df_mutation = exampleMutation,
#'                df_bam =  exampleBAM,
#'                df_mut_call = exampleMutCall,
#'                ref_genome = BSgenome.Hsapiens.UCSC.hg38::
#'                               BSgenome.Hsapiens.UCSC.hg38,
#'                sample_name = "H15-11943-1-T_TDv3",
#'                read_length = 151,
#'                adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
#'                adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
#'                short_homology_search_length = 4,
#'                progress_bar = "Y")
#' }
#' @export
fun_read_check <- function(df_mutation,
                          df_bam,
                          df_mut_call,
                          ref_genome,
                          sample_name,
                          read_length,
                          adapter_1,
                          adapter_2,
                          short_homology_search_length,
                          progress_bar) {
  Chr <- NULL
  Pos <- NULL
  
  if (length(df_mutation[, 1]) > 0) {
    # initialize
    msec <- NULL
    homology_search <- NULL
    mut_depth <- NULL
    length_flag <- 0
    chrom <- ""
    pre_search_length <- 0
    pre_search_length_default <- 0
    post_search_length_default <- 20
    post_search_length <- 20
    minimum_hairpin_length <- 15
    ref_width <- 150
    hairpin_search_length <- 50
    max_mutation_search <- 50
    neighbor_length <- 20
    laxness <- 1

    # analyze each somatic mutation
    for (i in seq_len(length(df_mutation[, 1]))) {
      pos_err <- 0
      indel_status <- 0
      indel_length <- 0
      with_indel_1 <- FALSE
      with_indel_2 <- FALSE
      mut_call <- logical(0)
      co_mut_pre <- 3
      co_mut_post <- 3
      penalty_pre <- 0
      penalty_post <- 0
      penalty_post <- 0
      near_indel_pre <- 0
      near_indel_pre_candidate <- 0
      near_indel_post <- 0
      near_indel_post_candidate <- 0
      mut_depth_tmp <- rep(0, 160)
      caution <- ""
      # extract mutation supporting reads
      if (df_mutation[i, "Chr"] != chrom) {
        chrom <- df_mutation[i, "Chr"]
        id_no <- df_bam$rname == chrom
        df_bam_qname <- df_bam$qname[id_no]
        df_bam_seq <- df_bam$seq[id_no]
        df_bam_strand <- df_bam$strand[id_no]
        df_bam_cigar <- df_bam$cigar[id_no]
        df_bam_qual <- df_bam$qual[id_no]
        df_bam_pos <- df_bam$pos[id_no]
      }
      mut_read <- df_mut_call %>%
        filter(Chr == df_mutation[i, "Chr"] & Pos == df_mutation[i, "Pos"])
      mut_type <- str_split(df_mutation[i, "Mut_type"], "-")[[1]][[2]]
      if (mut_type == "ins") {
        indel_length <- nchar(df_mutation[i, "Alt"]) - 1
        indel_status <- 1
        with_indel_1 <- TRUE
        with_indel_2 <- TRUE
      } else if (mut_type == "del") {
        indel_length <- nchar(df_mutation[i, "Ref"]) - 1
        indel_status <- 1
        with_indel_1 <- TRUE
        with_indel_2 <- TRUE
      }
      if (dim(mut_read)[[1]] == 0 & indel_status == 1) {
        for (tmp in 1:max_mutation_search) {
          if (dim(mut_read)[[1]] == 0) {
            mut_read <- df_mut_call %>%
              filter(Chr == df_mutation[i, "Chr"] &
                     Pos == (df_mutation[i, "Pos"] - tmp))
            pos_err <- tmp
          }
        }
      }
      if (dim(mut_read)[1] > 0) {
        mut_detail <- str_split(mut_read$Mut, pattern = ";")[[1]]
        mut_read_id_list <- str_split(mut_read$Mut_ID, pattern = ";")[[1]]
        if (indel_status == 0) {
          mut_call <- which(mut_detail == str_sub(df_mutation[i, "Alt"],
                                                 start = 1, end = 1))
        } else if (mut_type == "del") {
          mut_call <- which(mut_detail ==
                             str_replace(df_mutation[i, "Ref"],
                                         pattern = df_mutation[i, "Alt"],
                                         replacement = ".-"))
        } else if (mut_type == "ins") {
          mut_call <- which(mut_detail ==
                             str_replace(df_mutation[i, "Alt"],
                                         pattern = df_mutation[i, "Ref"],
                                         replacement = ".+"))
        }
      }
      # if mutation supporting reads exist
      if (length(mut_call) > 0) {
        pre_support_length <- 0
        post_support_length <- 0
        short_support_length <- 0
        pre_minimum_length <- read_length
        post_minimum_length <- read_length
        pre_rep_status <- 0
        post_rep_status <- 0
        homopolymer_status <- 0
        flag_hairpin <- 0
        hairpin_length <- 0
        soft_clipped_read <- 0
        total_read <- 0
        low_quality_base <- 0
        homology_search_tmp <- NULL
        distant_homology <- 0
        mutated_id <- mut_read_id_list[mut_call]
        mut_read_id <- str_sub(str_split(mutated_id, ",")[[1]],
                               start = 1,
                               end = -2)
        mut_read_strand <- str_sub(str_split(mutated_id, ",")[[1]],
                                  start = -1,
                                  end = -1)
        mut_read_strand <- gsub("r", "-", gsub("f", "+", mut_read_strand))
        ref_seq <- ref_genome[[df_mutation[i, "Chr"]]][
          (df_mutation[i, "Pos"] - pos_err - ref_width):
          (df_mutation[i, "Pos"] - pos_err + ref_width)]
        ref_indel <- c(ref_seq[1:ref_width], DNAString(df_mutation[i, "Alt"]),
                      ref_seq[(ref_width + nchar(df_mutation[i, "Ref"]) + 1):
                                (2 * ref_width + 1)])
        neighbor_seq <- df_mutation[i, "Neighborhood_sequence"]
        alt_length <- nchar(df_mutation[i, "Alt"])
        # sequence information around the mutation position
        setting <- fun_setting(pre = pre_search_length_default,
                              post = post_search_length_default,
                              neighbor_seq = neighbor_seq,
                              neighbor_length = neighbor_length,
                              alt_length = alt_length)
        pre_search_length <- setting[[1]]
        post_search_length <- setting[[2]]
        peri_seq_1 <- setting[[3]]
        peri_seq_2 <- setting[[4]]
        near_list_1 <- df_mut_call %>%
          filter(Chr == df_mutation[i, "Chr"] &
                 Pos >= (df_mutation[i, "Pos"] - pos_err - pre_search_length) &
                 Pos <= (df_mutation[i, "Pos"] + post_search_length) &
                 Pos != (df_mutation[i, "Pos"] - pos_err))
        near_list_2 <- df_mut_call %>%
          filter(Chr == df_mutation[i, "Chr"] &
                 Pos >= (df_mutation[i, "Pos"] - pos_err - post_search_length) &
                 Pos <= (df_mutation[i, "Pos"] + pre_search_length) &
                 Pos != (df_mutation[i, "Pos"] - pos_err))
        mut_near_1 <- dim(near_list_1)[1]
        mut_near_2 <- dim(near_list_2)[1]
        if (all(nchar(near_list_1$Ref) != nchar(near_list_1$Alt))) {
          with_indel_1 <- TRUE
        }
        if (all(nchar(near_list_2$Ref) != nchar(near_list_2$Alt))) {
          with_indel_2 <- TRUE
        }
        # short repeat around indel mutations
        if (mut_type == "ins") {
          rep_status <- fun_repeat_check(
            rep_a = DNAString(df_mutation[i, "Ref"]),
            rep_b = DNAString(df_mutation[i, "Alt"]),
            ref_seq = ref_seq,
            ref_width = ref_width,
            del = 0)
          pre_rep_status <- rep_status[[1]]
          post_rep_status <- rep_status[[2]]
          homopolymer_status <- rep_status[[3]]
        } else if (mut_type == "del") {
          rep_status <- fun_repeat_check(
            rep_a = DNAString(df_mutation[i, "Alt"]),
            rep_b = DNAString(df_mutation[i, "Ref"]),
            ref_seq = ref_seq,
            ref_width = ref_width,
            del = 1)
          pre_rep_status <- rep_status[[1]]
          post_rep_status <- rep_status[[2]]
          homopolymer_status <- rep_status[[3]]
        }
        for (depth in 1:(read_length + 1)) {
          mut_depth_tmp[depth] <-
            sum(df_bam_pos == (df_mutation[i, "Pos"] + 1 - depth))
        }
        if (progress_bar == "Y") {
          pb <- utils::txtProgressBar(min = 0,
                              max = max(1, length(mut_read_id)),
                              width = 20,
                              style = 3)
          pb_t <- ceiling(length(mut_read_id) / 100)
        }
        # analyze each mutation supporting read
        for (j in seq_len(length(mut_read_id))) {
          if (length_flag == 1) {
            setting <- fun_setting(pre = pre_search_length_default,
                                  post = post_search_length_default,
                                  neighbor_seq = neighbor_seq,
                                  neighbor_length = neighbor_length,
                                  alt_length = alt_length)
            pre_search_length <- setting[[1]]
            post_search_length <- setting[[2]]
            peri_seq_1 <- setting[[3]]
            peri_seq_2 <- setting[[4]]
          }
          length_flag <- 0
          mut_position <- 0
          flag_1 <- pre_search_length
          flag_2 <- pre_search_length
          # progress bar
          if (progress_bar == "Y") {
            utils::setTxtProgressBar(pb, j)
            if ((j - 1) %% pb_t == 0) {
              cat(paste(" Mutation:", i, "/",
                        length(df_mutation[, 1]), "      "))
            }
          }
          # specific read selection
          id_no <- df_bam_qname == mut_read_id[[j]]
          df_seq <- df_bam_seq[id_no]
          df_strand <- df_bam_strand[id_no]
          df_cigar <- df_bam_cigar[id_no]
          df_qual <- df_bam_qual[id_no]
          df_seq <- df_seq[df_strand == mut_read_strand[[j]]]
          df_cigar <- df_cigar[df_strand == mut_read_strand[[j]]]
          df_qual <- df_qual[df_strand == mut_read_strand[[j]]]
          if (length(df_seq) > 1) {
            df_cigar <- df_cigar[which.max(width(df_seq))]
            df_qual <- df_qual[which.max(width(df_seq)), ]
            df_seq <- df_seq[which.max(width(df_seq)), ]
          }
          df_seq <- df_seq[[1]]
          df_qual <- as.vector(asc(as.character(df_qual[1])))
          if (length(df_seq) > 0) {

            # determine mutation position in each read
            mutation_supporting_1 <-
              matchPattern(peri_seq_1,
                           df_seq,
                           max.mismatch = mut_near_1 + laxness,
                           min.mismatch = 0,
                           with.indels = with_indel_1,
                           fixed = FALSE)
            mutation_supporting_2 <-
              matchPattern(peri_seq_2,
                           df_seq,
                           max.mismatch = mut_near_2 + laxness,
                           min.mismatch = 0,
                           with.indels = with_indel_2,
                           fixed = FALSE)
            if (length(mutation_supporting_1) != 1 &
               length(mutation_supporting_2) != 1) {
              length_flag <- 1
              if (length(mutation_supporting_1) == 0 &
                 length(mutation_supporting_2) == 0) {
                for (Lax_1 in seq(3, 9, length = 3)) {
                  for (Lax_2 in 2:0) {
                    if (length(mutation_supporting_1) != 1 &
                       length(mutation_supporting_2) != 1) {
                      setting <- fun_setting(
                        pre = pre_search_length_default,
                        post = post_search_length_default - Lax_1 * laxness,
                        neighbor_seq = neighbor_seq,
                        neighbor_length = neighbor_length,
                        alt_length = alt_length)
                      pre_search_length <- setting[[1]]
                      post_search_length <- setting[[2]]
                      peri_seq_1 <- setting[[3]]
                      peri_seq_2 <- setting[[4]]
                      mutation_supporting_1 <- matchPattern(
                        peri_seq_1, df_seq,
                        max.mismatch = mut_near_1 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = with_indel_1,
                        fixed = FALSE)
                      mutation_supporting_2 <- matchPattern(
                        peri_seq_2,
                        df_seq,
                        max.mismatch = mut_near_2 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = with_indel_2,
                        fixed = FALSE)
                    }
                  }
                }
              }
              if (length(mutation_supporting_1) != 1 &
                 length(mutation_supporting_2) != 1) {
                for (Lax_1 in seq(5, 15, length = 3)) {
                  for (Lax_2 in 3:0) {
                    if (length(mutation_supporting_1) != 1 &
                       length(mutation_supporting_2) != 1) {
                      setting <- fun_setting(
                        pre = pre_search_length_default + Lax_1 * laxness,
                        post = post_search_length_default,
                        neighbor_seq = neighbor_seq,
                        neighbor_length = neighbor_length,
                        alt_length = alt_length)
                      pre_search_length <- setting[[1]]
                      post_search_length <- setting[[2]]
                      peri_seq_1 <- setting[[3]]
                      peri_seq_2 <- setting[[4]]
                      mutation_supporting_1 <- matchPattern(
                        peri_seq_1,
                        df_seq,
                        max.mismatch = mut_near_1 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = with_indel_1,
                        fixed = FALSE)
                      mutation_supporting_2 <- matchPattern(
                        peri_seq_2,
                        df_seq,
                        max.mismatch = mut_near_2 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = with_indel_2,
                        fixed = FALSE)
                    }
                  }
                }
              }
            }
            if (length(mutation_supporting_1) == 1) {
              mut_position <- min(
                length(df_seq),
                start(mutation_supporting_1) + pre_search_length)
              flag_1 <- hairpin_search_length
            } else if (length(mutation_supporting_2) == 1) {
              mut_position <- min(
                length(df_seq),
                start(mutation_supporting_2) + post_search_length)
              flag_2 <- hairpin_search_length
            }
            if (mut_position > 0) {
              total_read <- total_read + 1
              flag_hairpin_tmp <- 0
              check_hairpin <- 1

              # search co-mutations on neighbor
              if (mut_position > 10) {
                near_indel_pre_candidate <- near_indel_pre_candidate + 1
                co_mut_pre_tmp <- length(
                  matchPattern(
                    df_seq[(mut_position - 10):mut_position],
                    ref_indel[(ref_width - 9):(ref_width + 1)],
                    max.mismatch = 3,
                    min.mismatch = 0,
                    with.indels = FALSE,
                    fixed = TRUE))
                if (co_mut_pre_tmp == 0) {
                  near_indel_pre <- near_indel_pre + 1
                }
              }

              if (mut_position < (length(df_seq) - 10)) {
                near_indel_post_candidate <- near_indel_post_candidate + 1
                co_mut_post_tmp <- length(
                  matchPattern(
                    df_seq[mut_position:(mut_position + 10)],
                    ref_indel[(ref_width + 1):(ref_width + 11)],
                    max.mismatch = 3,
                    min.mismatch = 0,
                    with.indels = FALSE,
                    fixed = TRUE))
                if (co_mut_post_tmp == 0) {
                  near_indel_post <- near_indel_post + 1
                }
              }

              if (indel_status == 1) {
                comut_flag <- TRUE
                for (comut in 0:3) {
                  if (comut_flag == TRUE) {
                    co_mut_pre_tmp <- length(
                      matchPattern(df_seq[
                        max(1, (mut_position - indel_length - 1)):
                        min((mut_position + nchar(df_mutation[i, "Alt"]) + 4),
                              length(df_seq))],
                        ref_indel[(ref_width - indel_length):
                                (ref_width + nchar(df_mutation[i, "Alt"]) + 5)],
                         max.mismatch = comut,
                         min.mismatch = comut,
                         with.indels = FALSE,
                         fixed = TRUE))
                  }
                  if (co_mut_pre_tmp > 0) {
                    co_mut_pre <- min(co_mut_pre, comut)
                    comut_flag <- FALSE
                  }
                }
                comut_flag <- TRUE
                for (comut in 0:3) {
                  if (comut_flag == TRUE) {
                    co_mut_post_tmp <- length(
                      matchPattern(df_seq[
                        max(1, (mut_position - 4)):
                        min((mut_position + nchar(df_mutation[i, "Alt"]) +
                               indel_length - 1),
                            length(df_seq))],
                        ref_indel[(ref_width - 3):
                                  (ref_width + nchar(df_mutation[i, "Alt"]) +
                                     indel_length)],
                        max.mismatch = comut,
                        min.mismatch = comut,
                        with.indels = FALSE,
                        fixed = TRUE))
                  }
                  if (co_mut_post_tmp > 0) {
                    co_mut_post <- min(co_mut_post, comut)
                    comut_flag <- FALSE
                  }
                }
              } else{
                if ((mut_position - max(1, alt_length * 4 - 5)) > 0) {
                  comut_flag <- TRUE
                  for (comut in 0:3) {
                    if (comut_flag == TRUE) {
                      co_mut_pre_tmp <- length(
                        matchPattern(df_seq[
                          max(1, (mut_position - max(1, alt_length * 4 - 5))):
                            mut_position],
                          ref_indel[(ref_width -
                                       max(1, alt_length * 4 - 5) + 1):
                                      (ref_width + 1)],
                          max.mismatch = comut,
                          min.mismatch = comut,
                          with.indels = FALSE,
                          fixed = TRUE))
                    }
                    if (co_mut_pre_tmp > 0) {
                      co_mut_pre <- min(co_mut_pre, comut)
                      comut_flag <- FALSE
                    }
                  }
                  comut_flag <- TRUE
                  for (comut in 0:3) {
                    if (comut_flag == TRUE) {
                      co_mut_post_tmp <- length(
                        matchPattern(df_seq[
                          mut_position:
                            min(length(df_seq),
                                (mut_position + max(1, alt_length * 4 - 5)))],
                          ref_indel[(ref_width + 1):
                                      (ref_width +
                                         max(1, alt_length * 4 - 5) + 1)],
                          max.mismatch = comut,
                          min.mismatch = comut,
                          with.indels = FALSE,
                          fixed = TRUE))
                    }
                    if (co_mut_post_tmp > 0) {
                      co_mut_post <- min(co_mut_post, comut)
                      comut_flag <- FALSE
                    }
                  }
                }
              }

              # hairpin length calculation
              hairpin_seq <- reverseComplement(
                          df_seq[max(1, (mut_position - flag_1)):
                                 min(length(df_seq),
                                     (mut_position + flag_2 + alt_length - 1))])
              if (flag_1 == pre_search_length){
                if (minimum_hairpin_length <= length(hairpin_seq)) {
                  for (hair in minimum_hairpin_length:length(hairpin_seq)) {
                    if (check_hairpin == 1) {
                      hairpin_status <- fun_hairpin_check(
                                hairpin_seq[(length(hairpin_seq) - hair + 1):
                                            length(hairpin_seq)],
                                ref_seq[(ref_width + 1):(2 * ref_width + 1)],
                                hairpin_length,
                                hair)
                      hairpin_length <- hairpin_status[[1]]
                      check_hairpin <- hairpin_status[[2]]
                      flag_hairpin_tmp <- max(flag_hairpin_tmp, check_hairpin)
                    }
                  }
                }
              }
              if (flag_2 == pre_search_length){
                if (minimum_hairpin_length <= length(hairpin_seq)) {
                  for (hair in minimum_hairpin_length:length(hairpin_seq)) {
                    if (check_hairpin == 1) {
                      hairpin_status <- fun_hairpin_check(
                        hairpin_seq[1:hair],
                        ref_seq[1:(ref_width + 1)],
                        hairpin_length,
                        hair)
                      hairpin_length <- hairpin_status[[1]]
                      check_hairpin <- hairpin_status[[2]]
                      flag_hairpin_tmp <- max(flag_hairpin_tmp, check_hairpin)
                    }
                  }
                }
              }
              flag_hairpin <- flag_hairpin + flag_hairpin_tmp

              # supporting length calculation
              support_status <- fun_support(df_cigar,
                                     df_seq,
                                     mut_read_strand[[j]],
                                     adapter_1, adapter_2,
                                     mut_position,
                                     alt_length,
                                     indel_status)
              pre_support_length_tmp <- support_status[[1]]
              post_support_length_tmp <- support_status[[2]]
              soft_clipped_read_tmp <- support_status[[3]]

              # adjustment for consecutive snv
              if (mut_type == "snv") {
                post_support_length_tmp <-
                  post_support_length_tmp - alt_length + 1
              }

              # read quality check
              low_quality_base <-
                low_quality_base + sum(df_qual < 51) / length(df_seq)

              # save sequence for homologous region search
              if (flag_hairpin_tmp == 0 & flag_hairpin == 0) {
                if (indel_status == 1) {
                  pre_homology_search_seq <- df_seq[1:
                    min(length(df_seq),
                        mut_position + short_homology_search_length +
                                       post_rep_status + alt_length - 1)]
                  post_homology_search_seq <- df_seq[
                    max(1, mut_position - short_homology_search_length
                                        - pre_rep_status):
                    length(df_seq)]
                } else {
                  pre_homology_search_seq <- df_seq[1:
                    min(length(df_seq),
                     mut_position + short_homology_search_length + alt_length)]
                  post_homology_search_seq <- df_seq[
                    max(1, mut_position - short_homology_search_length):
                    length(df_seq)]
                }
                if (mut_read_strand[[j]] == "-") {
                  pre_homology_search_seq <- fun_hairpin_trimming(
                    pre_homology_search_seq,
                    mut_read_strand[[j]],
                    adapter_1, adapter_2)
                  pre_homology_search_seq <-
                    reverseComplement(pre_homology_search_seq)
                } else {
                  post_homology_search_seq <- fun_hairpin_trimming(
                    post_homology_search_seq,
                    mut_read_strand[[j]],
                    adapter_1, adapter_2)
                  post_homology_search_seq <-
                    reverseComplement(post_homology_search_seq)
                }
                homology_search_tmp <- rbind(homology_search_tmp,
                  data.frame(sample_name,
                             Chr = df_mutation[i, "Chr"],
                             Pos = df_mutation[i, "Pos"],
                             Ref = df_mutation[i, "Ref"],
                             Alt = df_mutation[i, "Alt"],
                             Direction = "pre",
                             Seq = as.character(pre_homology_search_seq)))
                homology_search_tmp <- rbind(homology_search_tmp,
                  data.frame(sample_name,
                             Chr = df_mutation[i, "Chr"],
                             Pos = df_mutation[i, "Pos"],
                             Ref = df_mutation[i, "Ref"],
                             Alt = df_mutation[i, "Alt"],
                             Direction = "post",
                             Seq = as.character(post_homology_search_seq)))
              }

              # summary
              soft_clipped_read <- soft_clipped_read + soft_clipped_read_tmp
              pre_support_length <- max(pre_support_length,
                                       pre_support_length_tmp)
              post_support_length <- max(post_support_length,
                                        post_support_length_tmp)
              short_support_length <-
                max(short_support_length,
                    min(pre_support_length_tmp, post_support_length_tmp))
              pre_minimum_length <- min(pre_minimum_length,
                                       pre_support_length_tmp)
              post_minimum_length <- min(post_minimum_length,
                                        post_support_length_tmp)
            }
          }
        }
        if (indel_status == 1) {
          penalty_pre <- 5 * co_mut_pre
          penalty_post <- 5 * co_mut_post
        } else {
          penalty_pre <- max(0, 4 * alt_length - 5) + 5 * co_mut_pre
          penalty_post <- max(0, 4 * alt_length - 5) + 5 * co_mut_post
        }

        if ((fun_zero(near_indel_pre, near_indel_pre_candidate) == 1 &
             pre_minimum_length >= 10) |
            (fun_zero(near_indel_post, near_indel_post_candidate) == 1 &
             post_minimum_length >= 10)) {
          pre_support_length <- min(pre_support_length + post_minimum_length,
                                    read_length)
          post_support_length <- min(post_support_length + pre_minimum_length,
                                     read_length)
          pre_minimum_length <- 0
          post_minimum_length <- 0
          if (indel_status == 1) {
            penalty_pre <- 0
            penalty_post <- 0
          } else {
            penalty_pre <- max(0, 4 * alt_length - 5)
            penalty_post <- max(0, 4 * alt_length - 5)
          }
          caution <- paste(caution,
                          "minimum_lengths are set to be 0",
                          "because of neighbor co-mutations,")
        }

        # data formatting
        msec_tmp <- df_mutation[i, ] %>% dplyr::mutate(
          read_length = read_length,
          mut_type = mut_type,
          alt_length = alt_length,
          total_read = total_read,
          soft_clipped_read = soft_clipped_read,
          flag_hairpin = flag_hairpin,
          hairpin_length = hairpin_length,
          pre_support_length = pre_support_length,
          post_support_length = post_support_length,
          short_support_length = short_support_length,
          pre_minimum_length = pre_minimum_length,
          post_minimum_length = post_minimum_length,
          low_quality_base_rate_under_q18 =
            fun_zero(low_quality_base, total_read),
          pre_rep_status = pre_rep_status,
          post_rep_status = post_rep_status,
          homopolymer_status = homopolymer_status,
          indel_status = indel_status,
          indel_length = indel_length,
          distant_homology = distant_homology,
          penalty_pre = penalty_pre,
          penalty_post = penalty_post,
          caution = caution
        )
        if (flag_hairpin == 0) {
          homology_search <- rbind(homology_search, homology_search_tmp)
        }
      } else{
        msec_tmp <- df_mutation[i, ] %>% dplyr::mutate(
          read_length = read_length,
          mut_type = mut_type,
          alt_length = 0,
          total_read = 0,
          soft_clipped_read = 0,
          flag_hairpin = 0,
          hairpin_length = 0,
          pre_support_length = 0,
          post_support_length = 0,
          short_support_length = 0,
          pre_minimum_length = 0,
          post_minimum_length = 0,
          low_quality_base_rate_under_q18 = 0,
          pre_rep_status = 0,
          post_rep_status = 0,
          homopolymer_status = 0,
          indel_status = 0,
          indel_length = 0,
          distant_homology = 0,
          penalty_pre = 0,
          penalty_post = 0,
          caution = ""
        )
      }
      mut_depth_tmp <- c(0, mut_depth_tmp)
      mut_depth_tmp <- t(cumsum(mut_depth_tmp))
      colnames(mut_depth_tmp) <- c("Zero", paste("Depth", 0:159, sep = ""))
      msec <- rbind(msec, msec_tmp)
      mut_depth <- rbind(mut_depth, mut_depth_tmp)
    }
    return(list(msec, homology_search, mut_depth))
  } else{
    return(list(df_mutation, NULL, NULL))
  }
}




# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
