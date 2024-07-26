#' Read check function.
#'
#' This function attempts to check the mutation profile in each read.
#'
#' @param df_mutation Mutation information.
#' @param df_bam Data from the BAM file.
#' @param ref_genome Reference genome for the data.
#' @param sample_name Sample name (character)
#' @param read_length The read length in the sequence.
#' @param adapter_1 The Read 1 adapter sequence of the library.
#' @param adapter_2 The Read 2 adapter sequence of the library.　
#' @param short_homology_search_length Small sequence for homology search.
#' @param min_homology_search Minimum length to define "homologous".
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
#' @importFrom Biostrings matchPattern
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings trimLRPatterns
#' @importFrom Biostrings width
#' @importFrom Biostrings start
#' @importFrom Biostrings end
#' @importFrom Biostrings nchar
#' @importFrom Biostrings subseq
#' @importFrom BiocGenerics as.data.frame
#' @examples
#' \dontrun{
#' data(exampleMutation)
#' data(exampleBam)
#' fun_read_check(df_mutation = exampleMutation,
#' df_bam = exampleBam,
#' ref_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#' sample_name = "sample",
#' read_length = 150,
#' adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
#' adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
#' short_homology_search_length = 4,
#' min_homology_search = 40,
#' progress_bar = "N"
#' )
#' }
#' @export
fun_read_check <- function(df_mutation,
                           df_bam,
                           ref_genome,
                           sample_name,
                           read_length,
                           adapter_1,
                           adapter_2,
                           short_homology_search_length,
                           min_homology_search,
                           progress_bar) {
  Chr <- NULL
  Pos <- NULL

  df_chr <- df_bam$rname
  if (!is.character(df_chr)) {
    df_chr <- as.character(df_chr)
  }
  if(str_sub(df_chr[[1]], start=1, end=3) != "chr"){
    GenomeInfoDb::seqlevelsStyle(df_chr) <- "UCSC"
  }
  if(str_sub(df_mutation$Chr[[1]], start=1, end=3) != "chr"){
    df_mutation$Chr <<- paste("chr", df_mutation$Chr, sep="")
  }

  if (length(df_mutation[, 1]) > 0) {
    # initialize
    msec <- NULL
    homology_search <- NULL
    mut_depth_pre <- NULL
    mut_depth_post <- NULL
    mut_depth_short <- NULL
    chrom <- ""
    pre_search_length <- 0
    pre_search_length_default <- 0
    post_search_length_default <- 20
    post_search_length <- 20
    minimum_hairpin_length <- 15
    ref_width <- 200
    hairpin_search_length_1 <- 5
    hairpin_search_length_2 <- 10
    neighbor_length <- 20
    laxness <- 1
    mut_type_all <- str_split(df_mutation[, "Mut_type"], "-", simplify = T)[,2]
    indel_status_all <- ifelse(mut_type_all == "snv", 0, 1)
    indel_length_all <- abs(nchar(df_mutation[, "Alt"]) - 
                              nchar(df_mutation[, "Ref"]))
    alt_all <- nchar(df_mutation[, "Alt"])
    # analyze each somatic mutation
    for (i in seq_len(length(df_mutation[, 1]))) {
      total_read <- 0
      pos_err <- 0
      indel_status <- 0
      indel_length <- 0
      mut_call <- 0
      mut_read <- matrix(logical(0))
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
      pre_support_length <- 0
      post_support_length <- 0
      short_support_length <- 0
      pre_farthest <- 0
      post_farthest <- 0
      pre_minimum_length <- 0
      post_minimum_length <- 0
      soft_clipped_read <- 0
      low_quality_base <- 0
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
        df_bam_isize_chr <- df_bam$isize[id_no]
        ref_genome_chr <- ref_genome[[chrom]]
        gc()
      }
      id_no <- (df_bam_pos_chr > (df_mutation[i, "Pos"] - 200) &
                df_bam_pos_chr < (df_mutation[i, "Pos"] + 1))
      df_bam_qname <- df_bam_qname_chr[id_no]
      df_bam_seq <- df_bam_seq_chr[id_no]
      df_bam_strand <- df_bam_strand_chr[id_no]
      df_bam_cigar <- df_bam_cigar_chr[id_no]
      df_bam_pos <- df_bam_pos_chr[id_no]
      df_bam_qual <- df_bam_qual_chr[id_no]
      df_bam_isize <- df_bam_isize_chr[id_no]
      mut_type <- mut_type_all[i]
      Mut_type <- ifelse(mut_type == "snv", "M",
                         ifelse(mut_type == "del", "D",
                                ifelse(mut_type == "ins", "I", "")))
      alt_length <- nchar(df_mutation[i, "Alt"])
      
      if (length(df_bam_pos) > 0) {
        cigar_num_all <- str_split(df_bam_cigar,
                                   "[:upper:]",
                                   simplify = T)
        cigar_num_all <- matrix(
          as.integer(cigar_num_all[,-dim(cigar_num_all)[2]]),
          ncol = dim(cigar_num_all)[2]-1)
        cigar_num_all[is.na(cigar_num_all)] = 0
        cigar_type_all <-  as.matrix(str_split(df_bam_cigar,
                                               "[:digit:]+",
                                               simplify = T))[,-1,drop=F]
        cigar_pos_all <- df_bam_pos
        cigar_seq_all <- df_bam_seq
        cigar_qname_all <- df_bam_qname
        cigar_isize_all <- df_bam_isize
        cigar_qual_all <- as.list(sapply(as.character(df_bam_qual),
               function(x) strtoi(charToRaw(x), 16L),
               simplify = FALSE))
        names(cigar_qual_all) <- NULL
        cigar_strand_all <- gsub("\\-", "r",
                             gsub("\\+", "f", df_bam_strand))
        cigar_read_pos <- matrix(1,
                                 nrow=dim(cigar_num_all)[1],
                                 ncol=dim(cigar_num_all)[2] + 1)
        cigar_genome_pos <- matrix(cigar_pos_all,
                                 nrow=dim(cigar_num_all)[1],
                                 ncol=dim(cigar_num_all)[2] + 1)
        
        for (k in seq_len(dim(cigar_type_all)[2])) {
          cigar_read_pos[,(k+1):(dim(cigar_type_all)[2] + 1)] =
            cigar_read_pos[,(k+1):(dim(cigar_type_all)[2] + 1)] +
            (cigar_type_all[,k] %in% c("M","I","S")) *
            cigar_num_all[,k]
          cigar_genome_pos[,(k+1):(dim(cigar_type_all)[2] + 1)] =
            cigar_genome_pos[,(k+1):(dim(cigar_type_all)[2] + 1)] +
            (cigar_type_all[,k] %in% c("M","D")) *
            cigar_num_all[,k]
        }
        indel_status <- indel_status_all[i]
        indel_length <- indel_length_all[i]
        alt <- alt_all[i] - indel_status
        indel_pos <- rowSums((Mut_type %in% c("D", "I")) * 
                            (cigar_read_pos - indel_status) *
                            (cigar_genome_pos ==
                              (df_mutation[i, "Pos"] + indel_status)) *
                            cbind((cigar_type_all == Mut_type), FALSE) *
                            cbind((cigar_num_all == indel_length), FALSE))
        indel_pos <- indel_pos *
          (as.character(
            subseq(cigar_seq_all,
                   ifelse(indel_pos == 0, 1, indel_pos),
                   ifelse(indel_pos +
                            (nchar(cigar_seq_all) > (indel_pos + alt)) *
                            alt == 0, 1, indel_pos +
                            (nchar(cigar_seq_all) > (indel_pos + alt)) *
                            alt))) == df_mutation[i, "Alt"])
        snv_pos <- rowSums(cbind((cigar_type_all == Mut_type), FALSE) *
                   ((cigar_genome_pos - df_mutation[i, "Pos"]) <= 0) *
                   cbind((((cigar_genome_pos[,-1,drop=F]) - df_mutation[i, "Pos"]) > 0),
                         FALSE) *
                   (cigar_read_pos -
                      cigar_genome_pos +
                      df_mutation[i, "Pos"]))
        snv_pos <- snv_pos *
                   (as.character(
                     subseq(cigar_seq_all,
                            ifelse(snv_pos == 0, 1, snv_pos),
                            ifelse((snv_pos +
                             (nchar(cigar_seq_all) > (snv_pos + alt - 1)) *
                             (alt - 1)) == 0, 1, (snv_pos +
                                (nchar(cigar_seq_all) > (snv_pos + alt - 1)) *
                                (alt - 1))))) ==
                   df_mutation[i, "Alt"])
        mut_pos <- indel_pos + snv_pos
        if (sum(mut_pos, na.rm = T) > 0){
          pre_supporting_length <- df_mutation[i, "Pos"] -
            cigar_genome_pos[,1] +
            indel_status
          post_supporting_length <- 
            cigar_genome_pos[,dim(cigar_genome_pos)[2]] -
            df_mutation[i, "Pos"] -
            (1 - indel_status) * alt
          
          mutated_seq <- (mut_pos > 0) &
            lapply((cigar_qual_all), mean)>=53
          if (sum(mutated_seq, na.rm = T) > 0) {
            mut_read_id <- cigar_qname_all[mutated_seq]
            mut_read_strand <- cigar_strand_all[mutated_seq]
            mut_position_cigar <- mut_pos[mutated_seq]
            mut_qual <- cigar_qual_all[mutated_seq]
            mut_seq <- cigar_seq_all[mutated_seq]
            mut_isize <- cigar_isize_all[mutated_seq]
            mut_cigar <- as.matrix(cigar_type_all[mutated_seq,])
            mut_qname <- cigar_qname_all[mutated_seq]
            mut_pre_supporting_length = pre_supporting_length[mutated_seq]
            mut_post_supporting_length = post_supporting_length[mutated_seq]
            mut_short_supporting_length =
              ifelse(mut_pre_supporting_length < mut_post_supporting_length,
                     mut_pre_supporting_length,
                     mut_post_supporting_length)
            
            covering_seq <- (pre_supporting_length >= 0 &
                               post_supporting_length >= 0)
            if (sum(covering_seq, na.rm = T) > 0) {
              pre_supporting_length <- pre_supporting_length[covering_seq]
              post_supporting_length <- post_supporting_length[covering_seq]
              short_supporting_length =
                ifelse(pre_supporting_length < post_supporting_length,
                       pre_supporting_length,
                       post_supporting_length)
    
              mut_call <- 1
              soft_clipped_read <- sum(rowSums(mut_cigar == "S") > 0, na.rm = T)
    
              for (depth in seq_len(length(pre_supporting_length))) {
                mut_depth_pre_tmp[min(201,pre_supporting_length[depth] + 2)] <-
                  mut_depth_pre_tmp[min(201,pre_supporting_length[depth] + 2)] + 1
                mut_depth_post_tmp[min(201,post_supporting_length[depth] + 2)] <-
                  mut_depth_post_tmp[min(201,post_supporting_length[depth] + 2)] + 1
                mut_depth_short_tmp[min(101,short_supporting_length[depth] + 2)] <-
                  mut_depth_short_tmp[min(101,short_supporting_length[depth] + 2)] + 1
              }
              
              pre_support_length <- max(mut_pre_supporting_length)
              post_support_length <- max(mut_post_supporting_length)
              short_support_length <- max(mut_short_supporting_length)
              pre_minimum_length <- min(mut_pre_supporting_length)
              post_minimum_length <- min(mut_post_supporting_length)
    
              post_farthest <- max(ifelse(
                mut_read_strand == "f" & mut_isize < 1000 & mut_isize > 0,
                  ifelse(mut_post_supporting_length > mut_isize - mut_pre_supporting_length - alt_length,
                         mut_post_supporting_length, mut_isize - mut_pre_supporting_length - alt_length),
                         mut_post_supporting_length))
              pre_farthest <- max(ifelse(
                mut_read_strand == "r" & mut_isize > -1000 & mut_isize < 0,
                  ifelse(mut_pre_supporting_length > - mut_isize - mut_post_supporting_length - alt_length,
                         mut_pre_supporting_length, - mut_isize - mut_post_supporting_length - alt_length),
                         mut_pre_supporting_length))
    
              # read quality check
              low_quality_base <- sum(unlist(lapply(mut_qual,
                                        function(x){sum(x<51)/length(x)})),na.rm = T)
              df_qual_pre <- lapply(mut_qual,
                                    function(x){x[max(1, mut_position_cigar - 10):
                                                  max(1, mut_position_cigar - 1)]})
              df_qual_post <- lapply(mut_qual,
                                     function(x){x[min(length(x), mut_position_cigar + 1):
                                                     min(length(x), mut_position_cigar + 10)]})
              pre_mutation_quality_score <- sum(unlist(df_qual_pre) < 51, na.rm = T)
              pre_mutation_quality_num <- length(unlist(df_qual_pre))
              post_mutation_quality_score <- sum(unlist(df_qual_post) < 51, na.rm = T)
              post_mutation_quality_num <- length(unlist(df_qual_post))
            }
          }
        }
      }
      # if mutation supporting reads exist
      if (mut_call == 1) {
        pre_rep_status <- 0
        post_rep_status <- 0
        pre_rep_short <- 0
        post_rep_short <- 0
        homopolymer_status <- 0
        flag_hairpin <- 0
        hairpin_length <- 0
        distant_homology <- 0
        search_status_1 <- 0
        search_status_2 <- 0
        indel_flag <- 0
        neighbor_seq <- df_mutation[i, "Neighborhood_sequence"]
        ref_seq <- ref_genome_chr[
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
        near_list_1 <- df_mutation %>%
          filter(Chr == df_mutation[i, "Chr"] &
                   Pos >= (df_mutation[i, "Pos"] - pos_err - pre_search_length) &
                   Pos <= (df_mutation[i, "Pos"] + post_search_length) &
                   Pos != (df_mutation[i, "Pos"] - pos_err))
        near_list_2 <- df_mutation %>%
          filter(Chr == df_mutation[i, "Chr"] &
                   Pos >= (df_mutation[i, "Pos"] - pos_err - post_search_length) &
                   Pos <= (df_mutation[i, "Pos"] + pre_search_length) &
                   Pos != (df_mutation[i, "Pos"] - pos_err))
        mut_near_1 <- dim(near_list_1)[[1]]
        mut_near_2 <- dim(near_list_2)[[1]]
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
        # save sequence for homologous region search
        if (indel_status == 1) {
          pre_homology_search_seq <-
            subseq(mut_seq,
                   1,ifelse(nchar(mut_seq) < min_homology_search, nchar(mut_seq),
                             ifelse(mut_position_cigar +
                                      short_homology_search_length +
                                      post_rep_short +
                                      alt_length - 1 < min_homology_search, min_homology_search,
                              ifelse(mut_position_cigar +
                                      short_homology_search_length +
                                      post_rep_short +
                                      alt_length - 1 > nchar(mut_seq),
                                     nchar(mut_seq),
                                     mut_position_cigar +
                                       short_homology_search_length +
                                       post_rep_short +
                                       alt_length - 1))))
          post_homology_search_seq <-
            subseq(mut_seq,
                    ifelse(nchar(mut_seq) < min_homology_search, 1,
                    ifelse(mut_position_cigar -
                             short_homology_search_length -
                             pre_rep_short < 1, 1,
                    ifelse(nchar(mut_seq) -
                             mut_position_cigar +
                             short_homology_search_length +
                             pre_rep_short < 39, nchar(mut_seq) - 39,
                           mut_position_cigar -
                             short_homology_search_length -
                             pre_rep_short))))
        } else if (indel_flag == 1) {
          pre_homology_search_seq <-
            subseq(mut_seq,
                   1,ifelse(nchar(mut_seq) < min_homology_search, nchar(mut_seq),
                            ifelse(mut_position_cigar +
                                     short_homology_search_length +
                                     post_rep_short +
                                     alt_length < min_homology_search, min_homology_search,
                                   ifelse(mut_position_cigar +
                                            short_homology_search_length +
                                            post_rep_short +
                                            alt_length > nchar(mut_seq),
                                          nchar(mut_seq),
                                          mut_position_cigar +
                                            short_homology_search_length +
                                            post_rep_short +
                                            alt_length))))
          post_homology_search_seq <-
            subseq(mut_seq,
                   ifelse(nchar(mut_seq) < min_homology_search, 1,
                          ifelse(mut_position_cigar -
                                   short_homology_search_length -
                                   pre_rep_short < 1, 1,
                                 ifelse(nchar(mut_seq) -
                                          mut_position_cigar +
                                          short_homology_search_length +
                                          pre_rep_short < 39, nchar(mut_seq) - 39,
                                        mut_position_cigar -
                                          short_homology_search_length -
                                          pre_rep_short))))
        } else {
         pre_homology_search_seq <-
            subseq(mut_seq,
                   1,ifelse(nchar(mut_seq) < min_homology_search, nchar(mut_seq),
                            ifelse(mut_position_cigar +
                                     short_homology_search_length +
                                     post_rep_short +
                                     alt_length < min_homology_search, min_homology_search,
                                   ifelse(mut_position_cigar +
                                            short_homology_search_length +
                                            post_rep_short +
                                            alt_length - 1 > nchar(mut_seq),
                                          nchar(mut_seq),
                                          mut_position_cigar +
                                            short_homology_search_length +
                                            post_rep_short +
                                            alt_length - 1))))
          post_homology_search_seq <-
            subseq(mut_seq,
                   ifelse(nchar(mut_seq) < min_homology_search, 1,
                     ifelse(mut_position_cigar -
                       short_homology_search_length < 1, 1,
                     ifelse(nchar(mut_seq) -
                              mut_position_cigar +
                        short_homology_search_length < 39, nchar(mut_seq) - 39,
                        mut_position_cigar -
                        short_homology_search_length))))
        }
        homology_search <- rbind(homology_search,
                data.frame(sample_name = sample_name,
                           Chr = df_mutation[i, "Chr"],
                           Pos = df_mutation[i, "Pos"],
                           Ref = df_mutation[i, "Ref"],
                           Alt = df_mutation[i, "Alt"],
                           Direction = "pre",
                           Seq = as.character(pre_homology_search_seq)),
                data.frame(sample_name = sample_name,
                           Chr = df_mutation[i, "Chr"],
                           Pos = df_mutation[i, "Pos"],
                           Ref = df_mutation[i, "Ref"],
                           Alt = df_mutation[i, "Alt"],
                           Direction = "post",
                           Seq = as.character(post_homology_search_seq)))
        

       if (progress_bar == "Y") {
          pb <- utils::txtProgressBar(min = 0,
                                      max = max(1, length(mut_read_id)),
                                      width = 20,
                                      style = 3)
          pb_t <- ceiling(length(mut_read_id) / 100)
        }
        # analyze each mutation supporting read
        for (j in seq_len(length(mut_read_id))) {
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
          df_seq <- mut_seq[[j]]
          df_strand <- mut_read_strand[j]
          df_cigar <- mut_cigar[j,]
          if (nchar(df_seq) > 20) {
            # determine mutation position in each read
            mut_position = mut_position_cigar[[j]]
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
              
            }
          }
        }
        if (indel_status == 1) {
          penalty_pre <- indel_length + 1 + pre_rep_status + 5 * co_mut_pre
          penalty_post <- indel_length + 1 + post_rep_status + 5 * co_mut_post
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
          pre_farthest = pre_farthest,
          post_farthest = post_farthest,
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
          pre_farthest = 0,
          post_farthest = 0,
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
      }
      mut_depth_pre_tmp <- t(cumsum(mut_depth_pre_tmp))
      mut_depth_post_tmp <- t(cumsum(mut_depth_post_tmp))
      mut_depth_short_tmp <- t(cumsum(mut_depth_short_tmp))
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
