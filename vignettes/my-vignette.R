## ----setting------------------------------------------------------------------
wd <- "~"
knitr::opts_chunk$set(collapse = TRUE,
                      fig.width = 12,
                      fig.height = 8,
                      echo = TRUE,
                      warning = FALSE,
                      message = TRUE,
                      comment = "#>")
options(rmarkdown.html_vignette.check_title = FALSE,
        show.error.messages = FALSE,
        warn = -1)

progress_bar <- "N"

## ----packages-----------------------------------------------------------------
library(MicroSEC)

