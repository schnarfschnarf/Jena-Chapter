fit_or_load <- function(
    name,
    fit_fun,
    models_dir = here::here("output", "models"),
    overwrite = FALSE,
    verbose = TRUE
) {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required.")
  }
  
  dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
  model_file <- file.path(models_dir, paste0(name, ".rds"))
  
  if (!overwrite && file.exists(model_file)) {
    if (isTRUE(verbose)) message("Loading model: ", name)
    return(readRDS(model_file))
  }
  
  if (isTRUE(verbose)) message("Fitting model: ", name)
  fit <- fit_fun()
  saveRDS(fit, model_file)
  fit
}