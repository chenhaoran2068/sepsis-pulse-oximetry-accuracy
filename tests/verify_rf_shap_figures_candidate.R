#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L,3L) ||
    (length(args)==3L && args[[3L]]!="bundle"))
  stop("Usage: RF_VALUES_MANIFEST.tsv DISPLAY_RUNS_DIRECTORY [bundle]", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
m <- fread(args[[1L]])
display_dir <- args[[2L]]
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (nrow(m) != 5L || !setequal(m$cohort, cohorts)) stop("RF_SHAP_QA_MANIFEST_INVALID", call. = FALSE)
checks <- 0L
for (cohort_value in cohorts) {
  dir <- m[cohort == cohort_value]$rf_dir
  figure_number <- match(cohort_value, cohorts) + 6L
  stem <- if(length(args)==3L) paste0("supplementary-figure-",figure_number) else
    paste0("display-sf",figure_number,"-sep19-r1")
  receipt_path <- file.path(display_dir, paste0(stem,"-receipt.tsv"))
  pdf_path <- file.path(display_dir, paste0(stem,".pdf"))
  if (!file.exists(pdf_path) || !file.exists(receipt_path)) stop("RF_SHAP_QA_DISPLAY_MISSING", call. = FALSE)
  receipt <- fread(receipt_path)
  summary <- fread(file.path(dir, "rf_shap_summary.tsv"))
  sample <- fread(file.path(dir, "rf_shap_sample.tsv"))
  values <- fread(file.path(dir, "rf_shap_feature_values.tsv"))
  grouped <- readRDS(file.path(dir, "rf_shap_grouped.rds"))
  if (nrow(receipt) != 1L || receipt$cohort[[1L]] != cohort_value ||
      receipt$figure_id[[1L]] != paste0("Supplementary Figure ", figure_number) ||
      receipt$shap_patient_n[[1L]] != nrow(sample) ||
      receipt$shap_feature_n[[1L]] != nrow(summary) ||
      receipt$shap_point_n[[1L]] != length(grouped) ||
      !identical(as.character(values$pair_uid), as.character(sample$pair_uid)) ||
      any(abs(colMeans(abs(grouped)) - summary$mean_abs_shap) > 1e-7)) {
    stop(paste("RF_SHAP_QA_SOURCE_DISPLAY_MISMATCH", cohort_value), call. = FALSE)
  }
  for (pair in list(c("PaO2", "pao2_candidate_value", "pao2_observed_n"),
                    c("pH", "ph_candidate_value", "ph_observed_n"),
                    c("PaCO2", "paco2_candidate_value", "paco2_observed_n"))) {
    n <- sum(is.finite(values[[pair[[2L]]]]))
    if (receipt[[pair[[3L]]]][[1L]] != n) {
      stop(paste("RF_SHAP_QA_DEPENDENCE_COUNT_MISMATCH", cohort_value, pair[[1L]]), call. = FALSE)
    }
    checks <- checks + 1L
  }
  checks <- checks + 6L
}
cat("RF_SHAP_FIGURES_INDEPENDENT_QA_PASS checks=", checks, "\n", sep = "")
