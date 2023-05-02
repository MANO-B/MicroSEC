#' Devide function without 0/0 errors
#'
#' This function attempts to devide without 0/0 errors.
#'
#' @param a,b Integers
#' @return a devided by b
fun_zero <- function(a, b)ifelse(b == 0, 0, a / b)
# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
