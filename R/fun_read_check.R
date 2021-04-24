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
#' @return list(msec, homology_search, mut_depth)
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom stringr str_split
#' @importFrom stringr str_sub
#' @importFrom stringr str_replace
#' @importFrom stringr str_count
#' @importFrom stringr str_dup
#' @importFrom BiocGenerics which
#' @importFrom BiocGenerics which.max
#' @importFrom BiocGenerics rbind
#' @importFrom BiocGenerics colnames
#' @importFrom gtools asc
#' @importFrom Biostrings matchPattern
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings trimLRPatterns
#' @importFrom Biostrings width
#' @importFrom Biostrings start
#' @importFrom Biostrings end
#' @importFrom Biostrings nchar
#' @importFrom BiocGenerics as.data.frame
#' @examples
#' \donttest{
#' fun_read_check(short_homology_search_length = 4)
#' }
#' @export
fun_read_check <- function(short_homology_search_length) {
  if (!exists("df_mutation")) {
    df_mutation <<- exampleMutation
  }
  if (!exists("df_bam")) {
    df_bam <<- exampleBAM
  }
  if (!exists("df_mut_call")) {
    df_mut_call <<- exampleMutCall
  }
  if (!exists("organism")) {
    organism <<- "hg38"
  }
  if (!exists("ref_genome")) {
    ref_genome <<- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
  if (!exists("sample_name")) {
    sample_name <<- "H15-11943-1-T_TDv3"
  }
  if (!exists("read_length")) {
    read_length <<- 151
  }
  if (!exists("adapter_1")) {
    adapter_1 <<- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  }
  if (!exists("adapter_2")) {
    adapter_2 <<- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  }
  if (!exists("progress_bar")) {
    progress_bar <<- "Y"
  }
  if (!exists("list_exist")) {
    list_exist <<- TRUE
  } else if (!list_exist) {
    df_mut_call <<- data.frame(logical(0))
  }
  Chr <- NULL
  Pos <- NULL
  options(show.error.messages = TRUE, warn = -1)
  #options(show.error.messages = FALSE, warn = -1)
  
  df_chr <- df_bam$rname
  if (!is.character(df_chr)) {
    df_chr <- as.character(df_chr)
  }
  GenomeInfoDb::seqlevelsStyle(df_chr) <- "UCSC"
  
  if (length(df_mutation[, 1]) > 0) {
    # initialize
    msec <- NULL
    homology_search <- NULL
    mut_depth_pre <- NULL
    mut_depth_post <- NULL
    mut_depth_short <- NULL
    length_flag <- 0
    chrom <- ""
    pre_search_length <- 0
    pre_search_length_default <- 0
    post_search_length_default <- 20
    post_search_length <- 20
    minimum_hairpin_length <- 15
    ref_width <- 200
    hairpin_search_length_1 <- 5
    hairpin_search_length_2 <- 10
    max_mutation_search <- 50
    neighbor_length <- 20
    laxness <- 1
    
    # analyze each somatic mutation
    for (i in seq_len(length(df_mutation[, 1]))) {
      total_read <- 0
      pos_err <- 0
      indel_status <- 0
      indel_length <- 0
      mut_call <- logical(0)
      mut_read <- matrix(logical(0))
      mut_read_id_list <- character(0)
      df_bam_pos <- logical(0)
      co_mut_pre <- 3
      co_mut_post <- 3
      penalty_pre <- 0
      penalty_post <- 0
      penalty_post <- 0
      near_indel_pre <- 0
      near_indel_pre_candidate <- 0
      near_indel_post <- 0
      near_indel_post_candidate <- 0
      pre_mutation_quality_score <- 0
      pre_mutation_quality_num <- 0
      post_mutation_quality_score <- 0
      post_mutation_quality_num <- 0
      mut_depth_pre_tmp <- rep(0, 201)
      mut_depth_post_tmp <- rep(0, 201)
      mut_depth_short_tmp <- rep(0, 101)
      caution <- ""
      # extract mutation supporting reads
      if (df_mutation[i, "Chr"] != chrom) {
        chrom <- df_mutation[i, "Chr"]
        id_no <- df_chr == chrom
        df_bam_qname_chr <- df_bam$qname[id_no]
        df_bam_seq_chr <- df_bam$seq[id_no]
        df_bam_strand_chr <- df_bam$strand[id_no]
        df_bam_cigar_chr <- df_bam$cigar[id_no]
        df_bam_qual_chr <- df_bam$qual[id_no]
        df_bam_pos_chr <- df_bam$pos[id_no]
        ref_genome <<- NULL
        gc()
        gc()
        fun_load_genome(organism)
      }
      mut_type <- str_split(df_mutation[i, "Mut_type"], "-")[[1]][[2]]
      if (list_exist) {
        if (mut_type == "snv") {
          mut_read <- df_mut_call %>%
            filter(Chr == df_mutation[i, "Chr"] & Pos == df_mutation[i, "Pos"])
          if (dim(mut_read)[1] > 0) {
            mut_detail <- str_split(mut_read$Mut, pattern = ";")[[1]]
            mut_read_id_list <- str_split(mut_read$Mut_ID, pattern = ";")[[1]]
            mut_call <- which(mut_detail == str_sub(df_mutation[i, "Alt"],
                                                    start = 1, end = 1))
          }
        } else if (mut_type == "ins") {
          indel_length <- nchar(df_mutation[i, "Alt"]) - 1
          indel_status <- 1
          for (tmp in 0:max_mutation_search) {
            if (length(mut_call) == 0) {
              mut_read <- df_mut_call %>%
                filter(Chr == df_mutation[i, "Chr"] &
                         Pos == (df_mutation[i, "Pos"] - tmp))
              pos_err <- tmp
              if (dim(mut_read)[[1]] > 0) {
                mut_detail <- str_split(mut_read$Mut, pattern = ";")[[1]]
                mut_read_id_list <- str_split(
                  mut_read$Mut_ID, pattern = ";")[[1]]
                mut_call <- which(mut_detail ==
                                    str_replace(df_mutation[i, "Alt"],
                                                pattern = df_mutation[i, "Ref"],
                                                replacement = ".+"))
              }
            }
          }
        } else if (mut_type == "del") {
          indel_length <- nchar(df_mutation[i, "Ref"]) - 1
          indel_status <- 1
          for (tmp in 0:max_mutation_search) {
            if (length(mut_call) == 0) {
              mut_read <- df_mut_call %>%
                filter(Chr == df_mutation[i, "Chr"] &
                         Pos == (df_mutation[i, "Pos"] - tmp))
              pos_err <- tmp
              if (dim(mut_read)[[1]] > 0) {
                mut_detail <- str_split(mut_read$Mut, pattern = ";")[[1]]
                mut_read_id_list <- str_split(
                  mut_read$Mut_ID, pattern = ";")[[1]]
                mut_call <- which(mut_detail ==
                                    str_replace(df_mutation[i, "Ref"],
                                                pattern = df_mutation[i, "Alt"],
                                                replacement = ".-"))
              }
            }
          }
        }
      } else {
        id_no <- (df_bam_pos_chr > (df_mutation[i, "Pos"] - 200) &
                    df_bam_pos_chr < (df_mutation[i, "Pos"] + 1))
        df_bam_qname <- df_bam_qname_chr[id_no]
        df_bam_seq <- df_bam_seq_chr[id_no]
        df_bam_strand <- df_bam_strand_chr[id_no]
        df_bam_cigar <- df_bam_cigar_chr[id_no]
        df_bam_pos <- df_bam_pos_chr[id_no]
        df_bam_qual <- df_bam_qual_chr[id_no]
        check_first <- TRUE
        mut_position_cigar <- logical(0)
        if (length(df_bam_pos) > 0) {
          if (mut_type == "snv") {
            alt <- nchar(df_mutation[i, "Alt"])
            for (depth in seq_len(length(df_bam_pos))) {
              if (df_bam_pos[depth] <= df_mutation[i, "Pos"]) {
                cigar_num <- as.integer(str_split(df_bam_cigar[depth],
                                                  "[:upper:]")[[1]])
                cigar_type <- str_split(df_bam_cigar[depth],
                                        "[:digit:]+")[[1]][-1]
                cigar_pos <- df_bam_pos[depth]
                cigar_seq <- df_bam_seq[depth][[1]]
                cigar_qname <- df_bam_qname[depth]
                cigar_strand <- gsub("\\-", "r",
                                     gsub("\\+", "f", df_bam_strand[depth]))
                cigar_qual <- as.vector(
                  asc(as.character(df_bam_qual[[depth]])))
                if (mean(cigar_qual) >= 53) {
                  read_pos <- 1
                  for (k in seq_len(length(cigar_type))) {
                    if ((cigar_pos <= df_mutation[i, "Pos"]) &
                        ((cigar_pos + cigar_num[k]) >=
                         df_mutation[i, "Pos"]) &
                        cigar_type[k] == "M") {
                      snv_pos <- read_pos + df_mutation[i, "Pos"] - cigar_pos
                      if (snv_pos > 0 &
                          (snv_pos + alt - 1) <= length(cigar_seq)) {
                        if (as.character(
                          cigar_seq[snv_pos:(snv_pos + alt - 1)]) ==
                          df_mutation[i, "Alt"]) {
                          if (check_first) {
                            mut_read_id_list <- paste(cigar_qname,
                                                      cigar_strand, sep = "")
                            check_first <- FALSE
                          } else {
                            mut_read_id_list <- paste(mut_read_id_list,
                                                      ",",
                                                      cigar_qname,
                                                      cigar_strand, sep ="")
                          }
                          mut_position_cigar <- c(mut_position_cigar, snv_pos) 
                          mut_call <- 1
                        }
                      }
                    }
                    if (cigar_type[k] == "D" | cigar_type[k] == "M") {
                      cigar_pos <- cigar_pos + cigar_num[k]
                    }
                    if (cigar_type[k] != "D" & cigar_type[k] != "H") {
                      read_pos <- read_pos + cigar_num[k]
                    }
                  }
                }
              }
            }
            mut_read_id_list = list(mut_read_id_list)
          } else if (mut_type == "ins") {
            indel_length <- nchar(df_mutation[i, "Alt"]) - 1
            indel_status <- 1
            alt <- nchar(df_mutation[i, "Alt"]) - 1
            mut_pos <- NULL
            for (depth in seq_len(length(df_bam_pos))) {
              if (df_bam_pos[depth] <= df_mutation[i, "Pos"]) {
                cigar_num <- as.integer(str_split(df_bam_cigar[depth],
                                                  "[:upper:]")[[1]])
                cigar_type <- str_split(df_bam_cigar[depth],
                                        "[:digit:]+")[[1]][-1]
                cigar_pos <- df_bam_pos[depth]
                cigar_qual <- as.vector(
                  asc(as.character(df_bam_qual[[depth]])))
                if (mean(cigar_qual) >= 53) {
                  read_pos <- 1
                  mut_pos_tmp <- NULL
                  for (k in seq_len(length(cigar_type))) {
                    if ((cigar_num[k] == alt) &
                        cigar_type[k] == "I") {
                      cigar_seq <- df_bam_seq[depth][[1]]
                      if (as.character(
                        cigar_seq[read_pos:(read_pos + alt - 1)]) ==
                        str_sub(df_mutation[i, "Alt"], start = 2, end = -1)) {
                        mut_pos_tmp <- c(mut_pos_tmp, cigar_pos)
                      }
                    }
                    if (cigar_type[k] == "D" | cigar_type[k] == "M") {
                      cigar_pos <- cigar_pos + cigar_num[k]
                    }
                    if (cigar_type[k] != "D" & cigar_type[k] != "H") {
                      read_pos <- read_pos + cigar_num[k]
                    }
                  }
                  if (cigar_pos >= (df_mutation[i, "Pos"] + alt + 1)) {
                    mut_pos <- c(mut_pos, mut_pos_tmp)
                  }
                }
              }
            }
            if (!is.null(mut_pos)){
              mut_pos <- as.integer(names(rev(sort(table(mut_pos))))[1])
              for (depth in seq_len(length(df_bam_pos))) {
                if (df_bam_pos[depth] <= df_mutation[i, "Pos"]) {
                  cigar_num <- as.integer(str_split(df_bam_cigar[depth],
                                                    "[:upper:]")[[1]])
                  cigar_type <- str_split(df_bam_cigar[depth],
                                          "[:digit:]+")[[1]][-1]
                  cigar_pos <- df_bam_pos[depth]
                  cigar_qual <- as.vector(
                    asc(as.character(df_bam_qual[[depth]])))
                  read_pos <- 1
                  cigar_qual <- as.vector(
                    asc(as.character(df_bam_qual[[depth]])))
                  if (mean(cigar_qual) >= 53) {
                    for (k in seq_len(length(cigar_type))) {
                      if ((cigar_pos == mut_pos) &
                          (cigar_num[k] == alt) &
                          cigar_type[k] == "I") {
                        cigar_seq <- df_bam_seq[depth][[1]]
                        cigar_qname <- df_bam_qname[depth]
                        cigar_strand <- gsub("\\-", "r",
                                             gsub("\\+", "f",
                                                  df_bam_strand[depth]))
                        if (as.character(
                          cigar_seq[read_pos:(read_pos + alt - 1)]) ==
                          str_sub(df_mutation[i, "Alt"],
                                  start = 2, end = -1)) {
                          if (check_first) {
                            mut_read_id_list <- paste(cigar_qname,
                                                      cigar_strand, sep = "")
                            check_first <- FALSE
                          } else {
                            mut_read_id_list <- paste(mut_read_id_list,
                                                      ",",
                                                      cigar_qname,
                                                      cigar_strand, sep ="")
                          }
                          mut_call <- 1
                          mut_position_cigar <- c(mut_position_cigar,
                                                  read_pos + df_mutation[i, "Pos"] - mut_pos) 
                        }
                      }
                      if (cigar_type[k] == "D" | cigar_type[k] == "M") {
                        cigar_pos <- cigar_pos + cigar_num[k]
                      }
                      if (cigar_type[k] != "D" & cigar_type[k] != "H") {
                        read_pos <- read_pos + cigar_num[k]
                      }
                    }
                  }
                }
              }
              mut_read_id_list = list(mut_read_id_list)
            }
          } else if (mut_type == "del") {
            indel_length <- nchar(df_mutation[i, "Ref"]) - 1
            indel_status <- 1
            alt <- nchar(df_mutation[i, "Ref"]) - 1
            mut_pos <- NULL
            for (depth in seq_len(length(df_bam_pos))) {
              if (df_bam_pos[depth] <= df_mutation[i, "Pos"]) {
                cigar_num <- as.integer(str_split(df_bam_cigar[depth],
                                                  "[:upper:]")[[1]])
                cigar_type <- str_split(df_bam_cigar[depth],
                                        "[:digit:]+")[[1]][-1]
                cigar_pos <- df_bam_pos[depth]
                cigar_qual <- as.vector(asc(as.character(df_bam_qual[[depth]])))
                mut_pos_tmp <- NULL
                if (mean(cigar_qual) >= 53) {
                  read_pos <- 1
                  for (k in seq_len(length(cigar_type))) {
                    if (cigar_num[k] == alt &
                        cigar_type[k] == "D") {
                      mut_pos_tmp <- c(mut_pos_tmp, cigar_pos)
                    }
                    if (cigar_type[k] == "D" | cigar_type[k] == "M") {
                      cigar_pos <- cigar_pos + cigar_num[k]
                    }
                    if (cigar_type[k] != "D" & cigar_type[k] != "H") {
                      read_pos <- read_pos + cigar_num[k]
                    }
                  }
                  if (cigar_pos >= (df_mutation[i, "Pos"] + alt + 1)) {
                    mut_pos <- c(mut_pos, mut_pos_tmp)
                  }
                }
              }
            }
            if (!is.null(mut_pos)){
              mut_pos <- as.integer(names(rev(sort(table(mut_pos))))[1])
              for (depth in seq_len(length(df_bam_pos))) {
                if (df_bam_pos[depth] <= df_mutation[i, "Pos"]) {
                  cigar_num <- as.integer(str_split(df_bam_cigar[depth],
                                                    "[:upper:]")[[1]])
                  cigar_type <- str_split(df_bam_cigar[depth],
                                          "[:digit:]+")[[1]][-1]
                  cigar_pos <- df_bam_pos[depth]
                  read_pos <- 1
                  cigar_qual <- as.vector(
                    asc(as.character(df_bam_qual[[depth]])))
                  if (mean(cigar_qual) >= 53) {
                    for (k in seq_len(length(cigar_type))) {
                      if (cigar_pos == mut_pos &
                          cigar_num[k] == alt &
                          cigar_type[k] == "D") {
                        cigar_qname <- df_bam_qname[depth]
                        cigar_strand <- gsub("\\-", "r",
                                             gsub("\\+", "f",
                                                  df_bam_strand[depth]))
                        if (check_first) {
                          mut_read_id_list <- paste(cigar_qname,
                                                    cigar_strand, sep = "")
                          check_first <- FALSE
                        } else {
                          mut_read_id_list <- paste(mut_read_id_list,
                                                    ",",
                                                    cigar_qname,
                                                    cigar_strand, sep ="")
                        }
                        mut_call <- 1
                        mut_position_cigar <- c(mut_position_cigar,
                                    read_pos + df_mutation[i, "Pos"] - mut_pos) 
                      }
                      if (cigar_type[k] == "D" | cigar_type[k] == "M") {
                        cigar_pos <- cigar_pos + cigar_num[k]
                      }
                      if (cigar_type[k] != "D" & cigar_type[k] != "H") {
                        read_pos <- read_pos + cigar_num[k]
                      }
                    }
                  }
                }
              }
              mut_read_id_list = list(mut_read_id_list)
            }
          }
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
        pre_rep_short <- 0
        post_rep_short <- 0
        homopolymer_status <- 0
        flag_hairpin <- 0
        hairpin_length <- 0
        soft_clipped_read <- 0
        low_quality_base <- 0
        homology_search_tmp <- NULL
        distant_homology <- 0
        search_status_1 <- 0
        search_status_2 <- 0
        indel_flag <- 0
        mutated_id <- mut_read_id_list[mut_call]
        mut_read_id <- str_sub(str_split(mutated_id, ",")[[1]],
                               start = 1,
                               end = -2)
        mut_read_strand <- str_sub(str_split(mutated_id, ",")[[1]],
                                   start = -1,
                                   end = -1)
        mut_read_strand <- gsub("r", "-", gsub("f", "+", mut_read_strand))
        id_no <- (df_bam_pos_chr > (df_mutation[i, "Pos"] - 200) &
                    df_bam_pos_chr < (df_mutation[i, "Pos"] + 1))
        df_bam_qname <- df_bam_qname_chr[id_no]
        df_bam_seq <- df_bam_seq_chr[id_no]
        df_bam_strand <- df_bam_strand_chr[id_no]
        df_bam_cigar <- df_bam_cigar_chr[id_no]
        df_bam_qual <- df_bam_qual_chr[id_no]
        df_bam_pos <- df_bam_pos_chr[id_no]
        neighbor_seq <- df_mutation[i, "Neighborhood_sequence"]
        alt_length <- nchar(df_mutation[i, "Alt"])
        ref_seq <- ref_genome[[df_mutation[i, "Chr"]]][
          (df_mutation[i, "Pos"] - ref_width):
            (df_mutation[i, "Pos"] + ref_width)]
        ref_indel <- c(ref_seq[1:ref_width], DNAString(df_mutation[i, "Alt"]),
                       ref_seq[(ref_width + nchar(df_mutation[i, "Ref"]) + 1):
                                 (2 * ref_width + 1)])
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
        # short repeat around indel mutations
        if (str_count(neighbor_seq, str_dup("A", 15)) > 0 |
            str_count(neighbor_seq, str_dup("T", 15)) > 0 |
            str_count(neighbor_seq, str_dup("G", 15)) > 0 |
            str_count(neighbor_seq, str_dup("C", 15)) > 0) {
          homopolymer_status <- 15
        }
        if (mut_type == "ins") {
          rep_status <- fun_repeat_check(
            rep_a = DNAString(df_mutation[i, "Ref"]),
            rep_b = DNAString(df_mutation[i, "Alt"]),
            ref_seq = ref_seq,
            ref_width = ref_width,
            del = 0)
          pre_rep_status <- rep_status[[1]]
          post_rep_status <- rep_status[[2]]
          pre_rep_short <- rep_status[[3]]
          post_rep_short <- rep_status[[4]]
          homopolymer_status <- max(homopolymer_status, rep_status[[5]])
        } else if (mut_type == "del") {
          rep_status <- fun_repeat_check(
            rep_a = DNAString(df_mutation[i, "Alt"]),
            rep_b = DNAString(df_mutation[i, "Ref"]),
            ref_seq = ref_seq,
            ref_width = ref_width,
            del = 1)
          pre_rep_status <- rep_status[[1]]
          post_rep_status <- rep_status[[2]]
          pre_rep_short <- rep_status[[3]]
          post_rep_short <- rep_status[[4]]
          homopolymer_status <- max(homopolymer_status, rep_status[[5]])
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
            df_qual <- df_qual[which.max(width(df_seq))]
            df_seq <- df_seq[which.max(width(df_seq))]
          }
          df_seq <- df_seq[[1]]
          df_qual <- as.vector(asc(as.character(df_qual[1])))
          if (length(df_seq) > 20) {
            # determine mutation position in each read
            if (list_exist) {
              mutation_supporting_1 <-
                matchPattern(peri_seq_1,
                             df_seq,
                             max.mismatch = mut_near_1,
                             min.mismatch = 0,
                             with.indels = FALSE,
                             fixed = FALSE)
              if (length(mutation_supporting_1) != 1) {
                length_flag <- 1
                for (Lax_1 in seq(0, 9, length = 4)) {
                  for (Lax_2 in 0:3) {
                    if (length(mutation_supporting_1) != 1) {
                      search_status_1 <- search_status_1 + 1
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
                        peri_seq_1,
                        df_seq,
                        max.mismatch = mut_near_1 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = FALSE,
                        fixed = FALSE)
                    }
                  }
                }
              }
              if (length(mutation_supporting_1) != 1) {
                for (Lax_1 in seq(5, 15, length = 3)) {
                  for (Lax_2 in 0:3) {
                    if (length(mutation_supporting_1) != 1) {
                      search_status_1 <- search_status_1 + 1
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
                        with.indels = FALSE,
                        fixed = FALSE)
                    }
                  }
                }
              }
              if (length(mutation_supporting_1) != 1) {
                length_flag <- 1
                if (length(mutation_supporting_1) == 0) {
                  for (Lax_1 in seq(0, 9, length = 4)) {
                    for (Lax_2 in 4:5) {
                      if (length(mutation_supporting_1) != 1) {
                        search_status_1 <- search_status_1 + 1
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
                          peri_seq_1,
                          df_seq,
                          max.mismatch = mut_near_1 + Lax_2 * laxness,
                          min.mismatch = 0,
                          with.indels = FALSE,
                          fixed = FALSE)
                      }
                    }
                  }
                }
                if (length(mutation_supporting_1) != 1) {
                  for (Lax_1 in seq(5, 20, length = 3)) {
                    for (Lax_2 in 4:5) {
                      if (length(mutation_supporting_1) != 1) {
                        search_status_1 <- search_status_1 + 1
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
                          with.indels = FALSE,
                          fixed = FALSE)
                      }
                    }
                  }
                }
              }
              if (length(mutation_supporting_1) != 1) {
                for (Lax_1 in seq(0, 9, length = 4)) {
                  for (Lax_2 in 0:3) {
                    if (length(mutation_supporting_1) != 1) {
                      search_status_1 <- search_status_1 + 1
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
                        peri_seq_1,
                        df_seq,
                        max.mismatch = mut_near_1 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = TRUE,
                        fixed = FALSE)
                    }
                  }
                }
              }
              pre_search_length_1 <- pre_search_length
              post_search_length_1 <- post_search_length
              pre_search_length <- pre_search_length_default
              post_search_length <- post_search_length_default
              mutation_supporting_2 <-
                matchPattern(peri_seq_2,
                             df_seq,
                             max.mismatch = mut_near_2,
                             min.mismatch = 0,
                             with.indels = FALSE,
                             fixed = FALSE)
              if (length(mutation_supporting_2) != 1) {
                length_flag <- 1
                for (Lax_1 in seq(0, 9, length = 4)) {
                  for (Lax_2 in 0:3) {
                    if (length(mutation_supporting_2) != 1) {
                      search_status_2 <- search_status_2 + 1
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
                      mutation_supporting_2 <- matchPattern(
                        peri_seq_2,
                        df_seq,
                        max.mismatch = mut_near_2 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = FALSE,
                        fixed = FALSE)
                    }
                  }
                }
              }
              if (length(mutation_supporting_2) != 1) {
                for (Lax_1 in seq(5, 15, length = 3)) {
                  for (Lax_2 in 0:3) {
                    if (length(mutation_supporting_2) != 1) {
                      search_status_2 <- search_status_2 + 1
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
                      mutation_supporting_2 <- matchPattern(
                        peri_seq_2,
                        df_seq,
                        max.mismatch = mut_near_2 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = FALSE,
                        fixed = FALSE)
                    }
                  }
                }
              }
              if (length(mutation_supporting_2) != 1) {
                for (Lax_1 in seq(0, 9, length = 4)) {
                  for (Lax_2 in 4:5) {
                    if (length(mutation_supporting_2) != 1) {
                      search_status_2 <- search_status_2 + 1
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
                      mutation_supporting_2 <- matchPattern(
                        peri_seq_2,
                        df_seq,
                        max.mismatch = mut_near_2 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = FALSE,
                        fixed = FALSE)
                    }
                  }
                }
              }
              if (length(mutation_supporting_2) != 1) {
                for (Lax_1 in seq(5, 20, length = 3)) {
                  for (Lax_2 in 4:5) {
                    if (length(mutation_supporting_2) != 1) {
                      search_status_2 <- search_status_2 + 1
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
                      mutation_supporting_2 <- matchPattern(
                        peri_seq_2,
                        df_seq,
                        max.mismatch = mut_near_2 + Lax_2 * laxness,
                        min.mismatch = 0,
                        with.indels = FALSE,
                        fixed = FALSE)
                    }
                  }
                }
              }
              if (length(mutation_supporting_1) == 1 &
                  length(mutation_supporting_2) == 1) {
                mut_position_1 <- min(
                  length(df_seq),
                  start(mutation_supporting_1) + pre_search_length_1)
                mut_position_2 <- min(
                  length(df_seq),
                  end(mutation_supporting_2) - pre_search_length - alt_length + 1)
                if (mut_position_1 == mut_position_2) {
                  mut_position <- mut_position_1
                } else if (mut_position_1 < mut_position_2) {
                  if (mut_position_1 > (mut_position_2 - 5) &
                      mut_position_1 < length(df_seq) &
                      str_count(as.character(as.data.frame(
                        mutation_supporting_1)[1]), "N") < 10 &
                      str_count(as.character(as.data.frame(
                        mutation_supporting_2)[1]), "N") < 10 &
                      search_status_1 < 4 & search_status_2 < 4) {
                    indel_flag <- 1
                    mut_position <- mut_position_1
                    rep_status <- fun_repeat_check(
                      df_seq[mut_position_1],
                      df_seq[mut_position_1:(mut_position_1 + 1)],
                      ref_seq,
                      ref_width - 1,
                      del = 1)
                    pre_rep_status <- max(pre_rep_status, rep_status[[1]])
                    post_rep_status <- max(post_rep_status, rep_status[[2]])
                    pre_rep_short <- max(pre_rep_short, rep_status[[3]])
                    post_rep_short <- max(post_rep_short, rep_status[[4]])
                    homopolymer_status <- max(homopolymer_status, rep_status[[3]])
                  } else if (str_count(as.character(as.data.frame(
                    mutation_supporting_1)[1]), "N") < 10 &
                    str_count(as.character(as.data.frame(
                      mutation_supporting_1)[1]), "N") < 10){
                    if (search_status_1 <= search_status_2) {
                      if (mut_position_1 > 0) {
                        mut_position <- mut_position_1
                      } else if (mut_position_2 <= length(df_seq)) {
                        mut_position <- mut_position_2
                      }
                    } else {
                      if (mut_position_2 <= length(df_seq)) {
                        mut_position <- mut_position_2
                      } else if (mut_position_1 > 0) {
                        mut_position <- mut_position_1
                      }
                    }
                  } else if (str_count(as.character(as.data.frame(
                    mutation_supporting_1)[1]), "N") < 10) {
                    if (mut_position_1 > 0) {
                      mut_position <- mut_position_1
                    }
                  } else if (str_count(as.character(as.data.frame(
                    mutation_supporting_2)[1]), "N") < 10) {
                    if (mut_position_2 <= length(df_seq)) {
                      mut_position <- mut_position_2
                    }
                  }
                } else if (mut_position_2 < mut_position_1 &
                           mut_position_2 < length(df_seq) &
                           mut_position_2 > (mut_position_1 - 5) &
                           str_count(as.character(as.data.frame(
                             mutation_supporting_1)[1]), "N") < 10 &
                           str_count(as.character(as.data.frame(
                             mutation_supporting_2)[1]), "N") < 10 &
                           search_status_1 < 4 & search_status_2 < 4) {
                  indel_flag <- 1
                  mut_position <- mut_position_2
                  rep_status <- fun_repeat_check(
                    df_seq[mut_position_2],
                    df_seq[mut_position_2:(mut_position_2 + 1)],
                    ref_seq,
                    ref_width,
                    del = 0)
                  pre_rep_status <- max(pre_rep_status, rep_status[[1]])
                  post_rep_status <- max(post_rep_status, rep_status[[2]])
                  pre_rep_short <- max(pre_rep_short, rep_status[[3]])
                  post_rep_short <- max(post_rep_short, rep_status[[4]])
                  homopolymer_status <- max(homopolymer_status, rep_status[[3]])
                } else if (str_count(as.character(as.data.frame(
                  mutation_supporting_1)[1]), "N") < 10 &
                  str_count(as.character(as.data.frame(
                    mutation_supporting_1)[1]), "N") < 10){
                  if (search_status_1 <= search_status_2) {
                    if (mut_position_1 > 0) {
                      mut_position <- mut_position_1
                    } else if (mut_position_2 <= length(df_seq)) {
                      mut_position <- mut_position_2
                    }
                  } else {
                    if (mut_position_2 <= length(df_seq)) {
                      mut_position <- mut_position_2
                    } else if (mut_position_1 > 0) {
                      mut_position <- mut_position_1
                    }
                  }
                } else if (str_count(as.character(as.data.frame(
                  mutation_supporting_1)[1]), "N") < 10) {
                  if (mut_position_1 > 0) {
                    mut_position <- mut_position_1
                  }
                } else if (str_count(as.character(as.data.frame(
                  mutation_supporting_2)[1]), "N") < 10) {
                  if (mut_position_2 <= length(df_seq)) {
                    mut_position <- mut_position_2
                  }
                }
              } else if (length(mutation_supporting_1) == 1) {
                if (str_count(as.character(as.data.frame(
                  mutation_supporting_1)[1]), "N") <= 10 &
                  start(mutation_supporting_1) > 0) {
                  mut_position <- min(
                    length(df_seq),
                    start(mutation_supporting_1) + pre_search_length_1)
                }
              } else if (length(mutation_supporting_2) == 1) {
                if (str_count(as.character(as.data.frame(
                  mutation_supporting_2)[1]), "N") <= 10 &
                  end(mutation_supporting_2) <= length(df_seq)) {
                  mut_position <- min(
                    length(df_seq),
                    end(mutation_supporting_2) -
                      pre_search_length - alt_length + 1)
                }
              }
            } else {
              mut_position = mut_position_cigar[[j]]
            }
            if (mut_position > 0 & mut_position <= length(df_seq)) {
              total_read <- total_read + 1
              flag_hairpin_tmp <- 0
              
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
              if (mut_position <= (length(df_seq) - 9 - alt_length)) {
                near_indel_post_candidate <- near_indel_post_candidate + 1
                co_mut_post_tmp <- length(
                  matchPattern(
                    df_seq[mut_position:(mut_position + 9 + alt_length)],
                    ref_indel[(ref_width + 1):(ref_width + 10 + alt_length)],
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
                df_seq[max(1, (mut_position - hairpin_search_length_1)):
                         min(length(df_seq),
                             (mut_position +
                                hairpin_search_length_2 +
                                alt_length - 1))])
              if (minimum_hairpin_length <= length(hairpin_seq)) {
                hairpin_status <- fun_hairpin_check(
                  hairpin_seq[(length(hairpin_seq) -
                                 minimum_hairpin_length + 1):
                                length(hairpin_seq)],
                  ref_seq,
                  hairpin_length,
                  minimum_hairpin_length)
                hairpin_length <- max(hairpin_length, hairpin_status[[1]])
                flag_hairpin_tmp <- max(flag_hairpin_tmp, hairpin_status[[2]])
              }
              hairpin_seq <- reverseComplement(
                df_seq[max(1, (mut_position - hairpin_search_length_2)):
                         min(length(df_seq),
                             (mut_position +
                                hairpin_search_length_1 +
                                alt_length - 1))])
              if (minimum_hairpin_length <= length(hairpin_seq)) {
                hairpin_status <- fun_hairpin_check(
                  hairpin_seq[1:minimum_hairpin_length],
                  ref_seq,
                  hairpin_length,
                  minimum_hairpin_length)
                hairpin_length <- max(hairpin_length, hairpin_status[[1]])
                flag_hairpin_tmp <- max(flag_hairpin_tmp, hairpin_status[[2]])
              }
              flag_hairpin <- flag_hairpin + flag_hairpin_tmp
              
              # supporting length calculation
              support_status <- fun_support(df_cigar,
                                            df_seq,
                                            mut_read_strand[[j]],
                                            adapter_1,
                                            adapter_2,
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
              df_qual_pre = df_qual[max(1, mut_position - 11):
                                      max(1, mut_position - 1)]
              df_qual_post = df_qual[min(length(df_seq), mut_position + 1):
                                       min(length(df_seq), mut_position + 10)]
              pre_mutation_quality_score <-
                pre_mutation_quality_score + sum(df_qual_pre < 51)
              pre_mutation_quality_num <-
                pre_mutation_quality_num + length(df_qual_pre)
              post_mutation_quality_score <-
                post_mutation_quality_score + sum(df_qual_post < 51)
              post_mutation_quality_num <-
                post_mutation_quality_num + length(df_qual_post)
              # save sequence for homologous region search
              if (flag_hairpin_tmp == 0 & flag_hairpin == 0) {
                if (indel_status == 1) {
                  pre_homology_search_seq <- df_seq[1:
                                                      min(length(df_seq),
                                                          mut_position + short_homology_search_length +
                                                            post_rep_short + alt_length - 1)]
                  post_homology_search_seq <- df_seq[
                    max(1, mut_position - short_homology_search_length
                        - pre_rep_short):
                      length(df_seq)]
                } else if (indel_flag == 1) {
                  pre_homology_search_seq <- df_seq[1:
                                                      min(length(df_seq),
                                                          mut_position + short_homology_search_length +
                                                            post_rep_short + alt_length)]
                  post_homology_search_seq <- df_seq[
                    max(1, mut_position - short_homology_search_length
                        - pre_rep_short):
                      length(df_seq)]
                } else {
                  pre_homology_search_seq <- df_seq[1:
                                                      min(length(df_seq),
                                                          mut_position + short_homology_search_length + alt_length)]
                  post_homology_search_seq <- df_seq[
                    max(1, mut_position - short_homology_search_length):
                      length(df_seq)]
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
        
        if (fun_zero(near_indel_pre, near_indel_pre_candidate) == 1 &
            pre_minimum_length >= 10){
          post_support_length <- min(post_support_length + pre_minimum_length,
                                     read_length)
          pre_minimum_length <- 0
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
        if (fun_zero(near_indel_post, near_indel_post_candidate) == 1 &
            post_minimum_length >= 10) {
          pre_support_length <- min(pre_support_length + post_minimum_length,
                                    read_length)
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
          low_quality_pre =
            fun_zero(pre_mutation_quality_score, pre_mutation_quality_num),
          low_quality_post =
            fun_zero(post_mutation_quality_score, post_mutation_quality_num),
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
      }
      if (total_read == 0) {
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
          low_quality_pre = 0,
          low_quality_post = 0,
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
      } else {
        # normal supporting status detection
        if (length(df_bam_pos) > 0) {
          for (depth in seq_len(length(df_bam_pos))) {
            if (df_bam_pos[depth] <= df_mutation[i, "Pos"]) {
              support_depth <- fun_normal_support(df_bam_cigar[depth],
                                                  df_bam_pos[depth],
                                                  df_mutation[i, "Pos"],
                                                  alt_length,
                                                  indel_status)
              if (support_depth[[3]] != -1) {
                mut_depth_pre_tmp[support_depth[[1]] + 2] <-
                  mut_depth_pre_tmp[support_depth[[1]] + 2] + 1
                mut_depth_post_tmp[support_depth[[2]] + 2] <-
                  mut_depth_post_tmp[support_depth[[2]] + 2] + 1
                mut_depth_short_tmp[support_depth[[3]] + 2] <-
                  mut_depth_short_tmp[support_depth[[3]] + 2] + 1
              }
            }
          }
        }
      }
      mut_depth_pre_tmp <- t(cumsum(mut_depth_pre_tmp))
      mut_depth_post_tmp <- t(cumsum(mut_depth_post_tmp))
      mut_depth_short_tmp <- t(cumsum(mut_depth_short_tmp))
      colnames(mut_depth_pre_tmp) <- c("Zero", paste("Depth", 0:199, sep = ""))
      colnames(mut_depth_post_tmp) <- c("Zero", paste("Depth", 0:199, sep = ""))
      colnames(mut_depth_short_tmp) <- c("Zero", paste("Depth", 0:99, sep = ""))
      mut_depth_pre <- rbind(mut_depth_pre, mut_depth_pre_tmp)
      mut_depth_post <- rbind(mut_depth_post, mut_depth_post_tmp)
      mut_depth_short <- rbind(mut_depth_short, mut_depth_short_tmp)
      msec <- rbind(msec, msec_tmp)
    }
    mut_depth <- list(mut_depth_pre, mut_depth_post, mut_depth_short)
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
