#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: RF_MANIFEST.tsv ROWS.tsv RECEIPT.tsv", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
m <- fread(args[[1L]])
rows <- fread(args[[2L]])
receipt <- fread(args[[3L]])
if (nrow(rows) != 5L || receipt$row_n[[1L]] != 5L) stop("RF_PERFORMANCE_ROW_COUNT_INVALID", call. = FALSE)
metrics <- c(rmse = "rmse_with_95ci", mae = "mae_with_95ci",
             r_squared = "r_squared_with_95ci",
             calibration_intercept = "calibration_intercept_with_95ci",
             calibration_slope = "calibration_slope_with_95ci")
checks <- 0L
for (i in seq_len(nrow(m))) {
  cohort_value <- m$cohort[[i]]
  row <- rows[cohort_key == cohort_value]
  dir <- m$rf_dir[[i]]
  source <- fread(file.path(dir, "rf_outer_test_performance.tsv"))
  hyper <- fread(file.path(dir, "rf_selected_hyperparameters.tsv"))
  if (nrow(row) != 1L || row$development_patients_pairs != paste0(source$development_patient_n[[1L]],
       "/", source$development_pair_n[[1L]]) ||
      row$heldout_patients_pairs != paste0(source$test_patient_n[[1L]], "/",
       source$test_pair_n[[1L]]) ||
      row$mtry_min_node_size != paste0(hyper$resolved_mtry_min[[1L]], "/",
       hyper$min_node_size[[1L]])) {
    stop(paste("RF_PERFORMANCE_COUNT_OR_HYPERPARAMETER_MISMATCH", cohort_value), call. = FALSE)
  }
  checks <- checks + 3L
  for (metric_value in names(metrics)) {
    z <- source[metric == metric_value]
    displayed <- sprintf("%.2f (%.2f to %.2f)", z$estimate, z$ci_lower, z$ci_upper)
    if (nrow(z) != 1L || row[[metrics[[metric_value]]]][[1L]] != displayed) {
      stop(paste("RF_PERFORMANCE_METRIC_MISMATCH", cohort_value, metric_value), call. = FALSE)
    }
    checks <- checks + 1L
  }
}
cat("RF_PERFORMANCE_TABLE_INDEPENDENT_QA_PASS checks=", checks, "\n", sep = "")
