#' Supporting length calculation function.
#'
#' This function attempts to calculate supporting lengths of a read.
#'
#' @param df_cigar The CIGAR data of the read.
#' @param df_seq The sequence to be checked.
#' @param mut_read_strand The strand of the sequence, "+" or "-".
#' @param ADAPTOR_SEQ The adapter sequence of the library.
#' @param mut_position The mutation position in the read.
#' @param Alt_length The length of altered bases.
#' @importFrom stringr str_split
#' @importFrom Biostrings trimLRPatterns
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @importFrom BiocGenerics as.data.frame
#' @return list(Pre_support_length, Post_support_length, Soft_Clipped_read)
#' @export
fun_support = function(df_cigar,
                       df_seq,
                       mut_read_strand,
                       ADAPTOR_SEQ,
                       mut_position,
                       Alt_length,
                       indel_status){
  Pre_support_length_tmp = 0
  Post_support_length_tmp = 0
  cigar_num = as.integer(str_split(df_cigar, "[:upper:]")[[1]])
  cigar_type = str_split(df_cigar, "[:digit:]+")[[1]][-1]
  tmp_pos = 1
  tmp_pos_indel = 0
  Soft_Clipped_read_tmp = 0
  for(k in 1:length(cigar_type)){
    if(cigar_type[k] == "H"){
      tmp_pos = tmp_pos + cigar_num[k]
    }
    if(cigar_type[k] == "I"){
      tmp_pos = tmp_pos + cigar_num[k]
    }
    if(cigar_type[k] == "S"){
      clipped_seq = df_seq[tmp_pos:(tmp_pos + cigar_num[k] - 1)]
      if(mut_read_strand == "+"){
        clipped_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ,
                                     subject = clipped_seq)
        clipped_seq = str_split(clipped_seq, ADAPTOR_SEQ)[[1]]
        if(length(clipped_seq)> 1){
          clipped_seq = clipped_seq[[1]]
        }
        clipped_seq = DNAString(clipped_seq)
      }
      else{
        clipped_seq = reverseComplement(clipped_seq)
        clipped_seq = trimLRPatterns(Rpattern = ADAPTOR_SEQ,
                                     subject = clipped_seq)
        clipped_seq = str_split(clipped_seq, ADAPTOR_SEQ)[[1]]
        if(length(clipped_seq)> 1){
          clipped_seq = clipped_seq[[1]]
        }
        clipped_seq = DNAString(clipped_seq)
      }
      if(length(clipped_seq) > 0){
        Soft_Clipped_read_tmp = 1
      }
      tmp_pos = tmp_pos + cigar_num[k]
    }
    if(cigar_type[k] == "M"){
      if((tmp_pos + cigar_num[k] - 1) < mut_position){
        Pre_support_length_tmp = 
          Pre_support_length_tmp + cigar_num[k]
      }
      else if(tmp_pos > (mut_position + Alt_length - 1)){
        Post_support_length_tmp = 
          Post_support_length_tmp + cigar_num[k]
      }
      else{
        if(indel_status == 1){
          Pre_support_length_tmp = 
            Pre_support_length_tmp + (mut_position - tmp_pos)
          Post_support_length_tmp = 
            Post_support_length_tmp + 
            (cigar_num[k] + tmp_pos - mut_position)
        }
        else{
          Pre_support_length_tmp =
            Pre_support_length_tmp + (mut_position - tmp_pos)
          Post_support_length_tmp = 
            Post_support_length_tmp + 
            (cigar_num[k] + tmp_pos - mut_position - Alt_length)
        }
      }
      tmp_pos = tmp_pos + cigar_num[k]
    }
  }
  return(list(Pre_support_length_tmp,
              Post_support_length_tmp,
              Soft_Clipped_read_tmp))
}

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
