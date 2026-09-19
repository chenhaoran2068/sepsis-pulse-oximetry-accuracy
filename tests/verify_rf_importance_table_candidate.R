#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: RF_MANIFEST.tsv ROWS.tsv RECEIPT.tsv", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
m <- fread(args[[1L]])
rows <- fread(args[[2L]])
receipt <- fread(args[[3L]])
if (nrow(rows) != 12L || receipt$feature_n[[1L]] != 12L) stop("RF_IMPORTANCE_OUTPUT_COUNT_INVALID", call. = FALSE)
checks <- 0L
total_sample <- 0L
for (i in seq_len(nrow(m))) {
  cohort_value <- m$cohort[[i]]
  shap <- fread(file.path(m$rf_dir[[i]], "rf_shap_summary.tsv"))
  sample <- fread(file.path(m$rf_dir[[i]], "rf_shap_sample.tsv"))
  total_sample <- total_sample + nrow(sample)
  for (j in seq_len(nrow(rows))) {
    feature_value <- rows$feature[[j]]
    z <- shap[feature_group == feature_value]
    if (nrow(z) == 0L) {
      if (rows[[cohort_value]][[j]] != "-") stop("MISSING_FEATURE_NOT_DASH", call. = FALSE)
    } else {
      higher <- sum(shap$mean_abs_shap > z$mean_abs_shap[[1L]])
      expected <- sprintf("%.2f (%d)", z$mean_abs_shap[[1L]], higher + 1L)
      if (rows[[cohort_value]][[j]] != expected) {
        stop(paste("RF_IMPORTANCE_VALUE_OR_RANK_MISMATCH", cohort_value, feature_value), call. = FALSE)
      }
    }
    checks <- checks + 1L
  }
}
if (total_sample != receipt$sample_n_total[[1L]]) stop("RF_IMPORTANCE_SAMPLE_TOTAL_MISMATCH", call. = FALSE)
cat("RF_IMPORTANCE_TABLE_INDEPENDENT_QA_PASS checks=", checks + 1L, "\n", sep = "")
