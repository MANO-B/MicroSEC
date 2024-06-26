% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_load_bam.R
\name{fun_load_bam}
\alias{fun_load_bam}
\title{BAM file loader}
\usage{
fun_load_bam(bam_file)
}
\arguments{
\item{bam_file}{Path of the BAM file.}
}
\value{
df_bam
}
\description{
This function attempts to load the BAM file.
}
\examples{
fun_load_bam(
  system.file("extdata", "sample.bam", package = "MicroSEC")
)
}
