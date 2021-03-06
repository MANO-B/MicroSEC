% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_homology.R
\name{fun_homology}
\alias{fun_homology}
\title{Homology check function.}
\usage{
fun_homology(msec, df_distant, min_homology_search)
}
\arguments{
\item{msec}{Mutation filtering information.}

\item{df_distant}{Sequences to be checked.}

\item{min_homology_search}{Minimum length to define "homologous".}
}
\value{
msec
}
\description{
This function attempts to search the homologous regions.
}
\examples{
\donttest{
fun_homology(msec = msec_read_checked,
             df_distant = homology_searched,
             min_homology_search = 40)
}
}
