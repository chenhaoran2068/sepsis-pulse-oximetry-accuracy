args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: COHORT ANALYSIS_INPUT_DIR FAILED_MI_OUTPUT_DIR", call. = FALSE)
cohort <- args[[1L]]
input <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mi <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)
extra_library <- Sys.getenv("R9_EXTRA_R_LIB", unset = "")
if (nzchar(extra_library)) .libPaths(c(extra_library, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(lme4); library(car)})
source(file.path(root, "code", "formal", "stage405_multilevel_core.R"))
source(file.path(root, "code", "formal", "stage405_formal_lmm_stream.R"))
x <- data.table::as.data.table(readRDS(file.path(input, "lmm_input_60.rds")))
coverage <- s405m_coverage_ledger(x, cohort)
predictors <- s405m_select_predictors(coverage)
d <- readRDS(file.path(mi, "imputations", "completed_01.rds"))
d <- s405l_prepare_model_data(d)
formula <- s405l_formula(cohort, predictors)
cat("FORMULA", paste(deparse(formula), collapse = " "), "\n")
for (optimizer in c("bobyqa", "nloptwrap")) {
  attempt <- s405l_fit_attempt(formula, d, optimizer)
  cat("OPTIMIZER", optimizer, "OK", attempt$ok,
      "MESSAGE", attempt$convergence_message, "\n")
  if (inherits(attempt$fit, "error")) cat("FIT_ERROR", conditionMessage(attempt$fit), "\n")
}
checkpoint_paths <- file.path(mi, "checkpoints", sprintf("fit_%02d.rds", 1:30))
if (all(file.exists(checkpoint_paths))) {
  fits <- lapply(checkpoint_paths, readRDS)
  pooled <- s405l_pool_stream_results(fits, data.table::uniqueN(d$cluster_patient))
  print(pooled$primary[, .(term, monte_carlo_error_to_se)])
}
