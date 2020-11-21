## -----------------------------------------------------------------------------
wd = "/mnt/HDD8TB/MicroSEC/archive"
knitr::opts_chunk$set(collapse = TRUE,
                      fig.width = 12,
                      fig.height = 8,
                      echo = TRUE,
                      warning = FALSE,
                      message = TRUE,
                      comment = "#>")
knitr::opts_knit$set(root.dir = wd)
options(rmarkdown.html_vignette.check_title = FALSE)
progress_bar = "N"

## ----packages-----------------------------------------------------------------
library(MicroSEC)

