#'MicroSEC: Microhomology-induced chimeric read-originating sequence error cleaning pipeline for FFPE samples
#'
#'MicroSEC sequence error filterinf pipeline for FFPE samples
#'Version 1.0.4: package in github (Oct. 31, 2020)
#'Version 1.0.3: parameter justification, version management with github (Oct. 29, 2020)
#'Version 1.0.2: Indel length correction based on Bowtie2 penalty score (Oct. 28, 2020)
#'Version 1.0.1: Correspond to the mouse genome. Homopolymer should have >=12 bases (Oct. 25, 2020)
#'Version 1.0.0: Original version (Oct. 23, 2020)
#'
#'Author: "Masachika Ikegami"
#'Date: "Oct. 28, 2020"
#'
#'This pipeline is designed for filtering mutations found in formalin-fixed and paraffin-embedded (FFPE) samples.
#'The MicroSEC filter utilizes a statistical analysis, and the results for mutations with less than 10 supporting reads are not reliable.
#'
#'Three files are nessesary for the analysis: mutation information file, BAM file, and mutation supporting read ID information file.
#'
#'File 1: mutation information file
#'This excel file should contain at least these contents:
#'      Sample     Gene HGVS.p Mut_type Total_QV>=20   %Alt  Chr       Pos Ref Alt SimpleRepeat_TRF                     Neighborhood_sequence  Transition
#'SL_1010-N6-B SLC25A24      _    1-snv          366 1.0929 chr1 108130741   C   T                N CTACCTGGAGAATGGGCCCATGTGTCCAGGTAGCAGTAAGC  C>T_t
#'Total_QV>=20: The read number with total Q-value >=20. 
#'SimpleRepeat_TRF: Whether the mutation locates at a simple repeat sequence or not.
#'Neighborhood_sequence: [5'-20 bases] + [Alt sequence] + [3'-20 bases]
#'Transition: 1-snv mutation pattern with a 3'-base. C>T_t represents CT to TT mutation. C>T_g_FFPE represents the possible FFPE artifact.
#'
#'File 2: BAM file
#'
#'File 3: mutation supporting read ID information file
#'This file should contain at least these contents:
#' Chr     Pos Ref Alt                                                                                                Mut_ID     Mut
#'chr1 2561609   T   A  _;ID001-1:579185f,ID004-1:1873933f;ID006-1:1131647f,ID001-1:570086f,ID008-1:1953407r,ID002-2:749570r  .;A;N#
#'
#'This pipeline contains 8 filtering processes.
#'
#'
#'Filter 1  : Shorter-supporting lengths distribute too short to occur (1-1 and 1-2).  
#'Filter 1-1: P-values are less than the threshold_p(default: 10^(-6)).  
#'Filter 1-2: The longest shorter-supporting lengths is shorter than 40% of the read length.  
#'Filter 2  : Palindromic sequences exist within 150 bases (2-1 or 2-2).  
#'Filter 2-1: A reverse complementary sequence consisting >= 30 bases exists on the opposite strand.  
#'Filter 2-2: >=50% mutation-supporting reads contains a reverse complementary sequence of the opposite strand consisting >= 15 bases.  
#'Filter 3  : 3’-/5’-supporting lengths are too densely distributed to occur (3-1, 3-2, and 3-3).  
#'Filter 3-1: P-values are less than the threshold_p(default: 10^(-6)).  
#'Filter 3-2: The distributions of 3’-/5’-supporting lengths are shorter than 80% of the read length.  
#'Filter 3-3: <10% of bases are low quality (Quality score <18).
#'Filter 4  : >=90% mutation-supporting reads are soft-clipped (after cutting adaptor sequence).  
#'Filter 5  : >=20% mutations were called by chimeric reads comprising two distant regions.
#'Filter 6  : Mutations locating at simple repeat sequences.
#'Filter 7  : C>T_g false positive calls in FFPE samples.
#'Filter 8  : Mutations locating at a >=15 homopolymer.
#'
#'Supporting lengths are adjusted considering small repeat sequences around the mutations.
#'
#'How to use in command line
#'Rscript MicroSEC.R [working/output directory] [sample name] [mutation information excel file] [BAM file] [read ID information directory] [read length] [adaptor sequence] [progress bar Y/N] [sample type: Human or Mouse]
#'
#'Example
#'Rscript MicroSEC.R /mnt/HDD8TB/FFPE SL_0002_L_FFPE_11-B /mnt/HDD8TB/FFPE/source/SL_mutation_somatic_200909.xlsx /mnt/HDD8TB/FFPE/source/SL_0002_L_FFPE_11_TDv4.realigned.bam /mnt/HDD8TB/FFPE/source/SL_mut_call_FFPE 150 AGATCGGAAGAGC Y Human
#'Rscript MicroSEC.R /mnt/HDD8TB/FFPE SL_1010-N6-B /mnt/HDD8TB/FFPE/source/SL_mutation_somatic_200907.xlsx /mnt/HDD8TB/FFPE/source/SL_1010-N6_TDv4tw.realigned.bam /mnt/HDD8TB/FFPE/source/SL_1010-N6-B 124 AGATCGGAAGAGC Y Human
#'Rscript MicroSEC.R /mnt/HDD8TB/FFPE PC9 /mnt/HDD8TB/FFPE/source/CCLE.xlsx /mnt/HDD8TB/FFPE/source/Cell_line/PC9_FEPE_Ag_TDv4.realigned.bam /mnt/HDD8TB/FFPE/source/PC9_FEPE 127 AGATCGGAAGAGC Y Human
#'
#'How to use in R console
#'
#'Sample name is set to the sample of interest in the Sample column of the mutation information file.
#'Confirm the read length in the plathome
#'Confirm the adaptor sequence; Todai Onco Panel ("AGATCGGAAGAGC")
#'If you want to know the progress visually, [progress bar Y/N] should be Y.
#'
#'Results are saved in a excel file.
#'The explatation of the results is written in detail in the second sheet of the excel file.




# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
