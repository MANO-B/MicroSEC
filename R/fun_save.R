#' Save function.
#'
#' This function attempts to save the filtering results.
#'
#' @param MSEC Mutation filtering information.
#' @param SAMPLE_INFO The sample name.
#' @param wd The directory to save.
#' @importFrom openxlsx createWorkbook
#' @importFrom openxlsx addWorksheet
#' @importFrom openxlsx writeData
#' @importFrom openxlsx freezePane
#' @importFrom openxlsx setColWidths
#' @importFrom openxlsx saveWorkbook
#' @export
fun_save = function(MSEC, SAMPLE_INFO, wd){
  # explanation
  MSEC_explain = data.frame(
    Total_read = "The number of reads supporting the mutation: [Total_read >= 10] filtering is strongly recommended",
    Soft_Clipped_read = "The number of soft-clipped reads",
    FLAG_Hairpin = "The number of reads containing a reverse complementary sequence of the opposite strand consisting >= 15 bases",
    Pre_support_length = "The longest mapped length of reads of the 5' sequence of the mutation",
    Post_support_length = "The longest mapped length of reads of the 3' sequence of the mutation",
    short_support_length = "The longest length of the shorter mapped length around the mutation of all mutation-supporting reads",
    distant_homology_rate = "The rate of chimeric reads comprising two distant regions",
    prob_Filter_1 = "The probability of the occurence of the shorter-supporting lengths distribution (threshold_p: 10^(-6))",
    prob_Filter_3_pre = "The probability of the occurence of the pre-supporting lengths distribution (threshold_p: 10^(-6))",
    prob_Filter_3_post = "The probability of the occurence of the post-supporting lengths distribution (threshold_p: 10^(-6))",
    Filter_1_mutation_intra_hairpin_loop = "Shorter-supporting lengths are too short to occur",
    Filter_2_hairpin_structure = "Palindromic sequences exist within 150 bases",
    Filter_3_microhomology_induced_mutation = "3’-/5’-supporting lengths are too densely distributed to occur.",
    Filter_4_soft_clipping = ">=90% mutation-supporting reads are soft-clipped",
    Filter_5_highly_homologous_region = ">=20% mutations were called by chimeric reads comprising two distant regions",
    Filter_6_simple_repeat = "Mutations locating at simple repeat sequences",
    Filter_7_C_to_G_artifact = "C>T_g false positive calls in FFPE samples",
    Filter_8_mutation_at_homopolymer =  "Indel mutations located inside a >=15 homopolymer",
    MSEC_filter_1234 = "One or more filters of Fitler 1, 2, 3, or 4 are TRUE",
    MSEC_filter_12345 = "One or more filters of Fitler 1, 2, 3, 4, or 5 are TRUE",
    MSEC_filter_all= "One or more filters are TRUE"
  )
  MSEC_explain = data.frame(t(MSEC_explain))
  colnames(MSEC_explain) = c("Explanation")
  MSEC_explain_Name = data.frame(rownames(MSEC_explain))
  colnames(MSEC_explain_Name) = c("Name")
  MSEC_explain = cbind(MSEC_explain_Name, MSEC_explain)
  colnames(MSEC_explain) = c("Name", "Explanation")
  
  # save the results
  NewWb = createWorkbook()
  addWorksheet(wb = NewWb, sheetName = "MicroSEC_results", gridLines = TRUE)
  addWorksheet(wb = NewWb, sheetName = "MicroSEC_explanation", gridLines = TRUE)
  if(!is.null(MSEC)){
    writeData(wb = NewWb, sheet = "MicroSEC_results", x = MSEC,
              xy = c(1,1), borders = "all", withFilter=TRUE)
  } else{
    writeData(wb = NewWb, sheet = "MicroSEC_results", x = "",
              xy = c(1,1), borders = "all", withFilter=TRUE)
  }
  writeData(wb = NewWb, sheet = "MicroSEC_explanation", x = MSEC_explain,
            xy = c(1,1), borders = "all")
  freezePane(wb = NewWb, sheet = "MicroSEC_results", firstActiveRow = 2)
  if(!is.null(MSEC)){
    setColWidths(wb = NewWb, sheet = "MicroSEC_results", cols = 1:ncol(MSEC), widths = "auto")
  }
  setColWidths(wb = NewWb, sheet = "MicroSEC_explanation", cols = 1:ncol(MSEC_explain), widths = "auto")
  
  saveWorkbook(wb = NewWb, file = paste(wd, "/MicroSEC-result_", SAMPLE_INFO, "_", Sys.Date(), ".xlsx", sep=""), overwrite = TRUE, returnValue = FALSE)
  if(!is.null(MSEC)){
    write.table(MSEC, file=paste(wd, "/MicroSEC_", SAMPLE_INFO, ".tsv", sep=""), sep = "\t", na="", row.names=FALSE, col.names=FALSE, quote=FALSE)
  } else{
    write.table("", file=paste(wd, "/MicroSEC_", SAMPLE_INFO, "_NULL.txt", sep=""), sep = "\t", na="", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }
}


# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
