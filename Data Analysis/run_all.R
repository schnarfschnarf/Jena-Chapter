#!/usr/bin/env Rscript
if (requireNamespace("renv", quietly = TRUE)) {
  renv::activate()
}

# Default: do NOT refit models
REFIT_MODELS <- identical(Sys.getenv("REFIT"), "1")

set.seed(1234)
options(mc.cores = 1)

library(here)
library(rmarkdown)

dir.create(here("3_Graphs"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "models"), showWarnings = FALSE, recursive = TRUE)

message("REFIT_MODELS = ", REFIT_MODELS)

rmarkdown::render(here("2_Code", "Data_Analysis_Economic.Rmd"), clean = TRUE)
rmarkdown::render(here("2_Code", "Data_Analysis_Ecological.Rmd"), clean = TRUE)
rmarkdown::render(here("2_Code", "Insurance_Analysis.Rmd"), clean = TRUE)
rmarkdown::render(here("2_Code", "Data_Analysis_Dominance.Rmd"), clean = TRUE)