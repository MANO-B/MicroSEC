#' Save function.
#'
#' This function attempts to save the filtering results.
#'
#' @param msec Mutation filtering information.
#' @param sample_info The sample name.
#' @param wd The directory to save.
#' @importFrom openxlsx createWorkbook
#' @importFrom openxlsx addWorksheet
#' @importFrom openxlsx writeData
#' @importFrom openxlsx freezePane
#' @importFrom openxlsx setColWidths
#' @importFrom openxlsx saveWorkbook
#' @examples
#' \donttest{
#' fun_save(msec_analyzed, "test_data", ".")
#' }
#' @export
fun_save <- function(msec, sample_info, wd) {
  # explanation
  msec_explain <- data.frame(
    total_read = paste(
      "The number of reads supporting the mutation:",
      "[Total_read >= 10] filtering is strongly recommended"),
    soft_clipped_read = "The number of soft-clipped reads",
    flag_hairpin = paste(
      "The number of reads containing a reverse complementary sequence of the",
      "opposite strand consisting >= 15 bases"),
    pre_support_length =
      "The longest mapped length of reads of the 5' sequence of the mutation",
    post_support_length =
      "The longest mapped length of reads of the 3' sequence of the mutation",
    short_support_length = paste(
      "The longest length of the shorter mapped length around the mutation of",
      "all mutation-supporting reads"),
    distant_homology_rate =
      "The rate of chimeric reads comprising two distant regions",
    low_quality_base_rate_under_q18 =
      "The rate of low quality bases in the mutation supporting reads",
    prob_filter_1 = paste(
      "The probability of the occurence of the shorter-supporting lengths",
      "distribution (threshold_p: 10 ^ (-6))"),
    prob_filter_3_pre = paste(
      "The probability of the occurence of the pre-supporting lengths",
      "distribution (threshold_p: 10 ^ (-6))"),
    prob_filter_3_post = paste(
      "The probability of the occurence of the post-supporting lengths",
      "distribution (threshold_p: 10 ^ (-6))"),
    filter_1_mutation_intra_hairpin_loop =
      "Shorter-supporting lengths are too short to occur",
    filter_2_hairpin_structure =
      "Palindromic sequences exist within 150 bases",
    filter_3_microhomology_induced_mutation =
      "3'-/5'-supporting lengths are too densely distributed to occur.",
    filter_4_soft_clipping = ">=90% mutation-supporting reads are soft-clipped",
    filter_5_highly_homologous_region = paste(
      ">=20% mutations were called by chimeric reads comprising two distant",
      "regions"),
    filter_6_simple_repeat = "Mutations locating at simple repeat sequences",
    filter_7_c_to_t_artifact = "C>T_g false positive calls in FFPE samples",
    filter_8_mutation_at_homopolymer =
      "Indel mutations located inside a >=15 homopolymer",
    msec_filter_1234 = "One or more filters of Fitler 1, 2, 3, or 4 are TRUE",
    msec_filter_12345 =
      "One or more filters of Fitler 1, 2, 3, 4, or 5 are TRUE",
    msec_filter_all = "One or more filters are TRUE"
  )
  msec_explain <- data.frame(t(msec_explain))
  colnames(msec_explain) <- c("Explanation")
  msec_explain_name <- data.frame(rownames(msec_explain))
  colnames(msec_explain_name) <- c("Name")
  msec_explain <- cbind(msec_explain_name, msec_explain)
  colnames(msec_explain) <- c("Name", "Explanation")

  # save the results
  new_wb <- createWorkbook()
  addWorksheet(wb = new_wb, sheetName = "MicroSEC_results",
               gridLines = TRUE)
  addWorksheet(wb = new_wb, sheetName = "MicroSEC_explanation",
               gridLines = TRUE)
  writeData(wb = new_wb, sheet = "MicroSEC_results", x = msec,
            xy = c(1, 1), borders = "all", withFilter = TRUE)
  writeData(wb = new_wb, sheet = "MicroSEC_explanation", x = msec_explain,
            xy = c(1, 1), borders = "all")
  freezePane(wb = new_wb, sheet = "MicroSEC_results", firstActiveRow = 2)
  setColWidths(wb = new_wb, sheet = "MicroSEC_results",
               cols = seq_len(ncol(msec)), widths = "auto")
  setColWidths(wb = new_wb, sheet = "MicroSEC_explanation",
               cols = seq_len(ncol(msec_explain)), widths = "auto")
  dir.create(file.path(wd), showWarnings = FALSE)
  saveWorkbook(wb = new_wb,
               file = paste(wd, "/MicroSEC-result_", sample_info, "_",
                            Sys.Date(), ".xlsx", sep = ""),
               overwrite = TRUE, returnValue = FALSE)
  utils::write.table(msec,
              file = paste(wd, "/MicroSEC_", sample_info, ".tsv", sep = ""),
              sep = "\t", na = "", row.names = FALSE, col.names = TRUE,
              quote = FALSE)
}


# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
