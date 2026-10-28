suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
qc_dir <- file.path(root, "results", "qc")
super_pop_shift <- c(AFR = -1, AMR = 0, EAS = 1, EUR = 0.5, SAS = -0.5)  # scenario B environment
VAR_G <- 0.4   # genetic variance share, both scenarios
VAR_S <- 0.2   # environmental (population-linked) share, scenario B only

lambda_gc <- function(p) median(qchisq(p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
