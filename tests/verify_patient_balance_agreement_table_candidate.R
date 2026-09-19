#!/usr/bin/env Rscript
# Reuse isolated invented-cohort agreement outputs, never approved clinical results.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: BA_MANIFEST60.tsv CANDIDATE_ROOT RENDERER.R NEW_TEST_DIR", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
source_manifest <- fread(args[[1L]])
candidate <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
renderer <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
test_dir <- args[[4L]]
if (file.exists(test_dir)) stop("TEST_DIR_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
dir.create(test_dir, recursive = TRUE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
m <- rbindlist(lapply(cohorts, function(cohort_value) {
  original <- source_manifest[cohort == cohort_value]
  if (nrow(original) != 1L) stop("SYNTHETIC_SOURCE_MANIFEST_INVALID", call. = FALSE)
  path <- file.path(candidate, "runs", paste0("agreement-extended-contract-",
      tolower(cohort_value), "60-sep19-r1"), "agreement_overall.tsv")
  data.table(cohort = cohort_value, analysis_rds = original$analysis_rds,
             agreement_tsv = path)
}))
manifest_path <- file.path(test_dir, "manifest.tsv")
fwrite(m, manifest_path, sep = "\t")
pdf_path <- file.path(test_dir, "table.pdf")
rows_path <- file.path(test_dir, "rows.tsv")
receipt_path <- file.path(test_dir, "receipt.tsv")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
out <- suppressWarnings(system2(rscript,
  c(shQuote(renderer), shQuote(manifest_path), shQuote(pdf_path), shQuote(rows_path),
    shQuote(receipt_path)), stdout = TRUE, stderr = TRUE))
status <- attr(out, "status")
if ((!is.null(status) && status != 0L) || !all(file.exists(c(pdf_path, rows_path, receipt_path)))) {
  stop(paste("BALANCE_TABLE_RENDER_FAILED", paste(out, collapse = " | ")), call. = FALSE)
}
rows <- fread(rows_path)
if (nrow(rows) != 5L) stop("BALANCE_TABLE_ROW_COUNT_INVALID", call. = FALSE)
fmt <- function(z, name) sprintf("%.2f (%.2f to %.2f)", z[[name]],
  z[[paste0(name, "_ci_lower")]], z[[paste0(name, "_ci_upper")]])
checks <- 0L
for (i in seq_along(cohorts)) {
  cohort_value <- cohorts[[i]]
  input_row <- m[cohort == cohort_value]
  o <- fread(input_row$agreement_tsv)
  b <- fread(file.path(dirname(input_row$agreement_tsv), "agreement_patient_balanced.tsv"))
  shown <- rows[i]
  if (shown$patients_n != format(o$patient_n, big.mark = ",", trim = TRUE) ||
      shown$balanced_mean_bias_with_95ci != fmt(b, "balanced_mean_bias") ||
      shown$balanced_lower_loa_with_95ci != fmt(b, "balanced_lower_loa") ||
      shown$balanced_upper_loa_with_95ci != fmt(b, "balanced_upper_loa") ||
      shown$balanced_arms_with_95ci != fmt(b, "balanced_arms") ||
      shown$empirical_p2_5_with_95ci != fmt(o, "empirical_p2_5") ||
      shown$empirical_p97_5_with_95ci != fmt(o, "empirical_p97_5") ||
      shown$patient_variance != sprintf("%.2f", o$between_patient_variance) ||
      shown$residual_variance != sprintf("%.2f", o$within_stay_residual_variance) ||
      shown$total_difference_sd != sprintf("%.2f", o$total_difference_sd)) {
    stop(paste("BALANCE_TABLE_DISPLAY_MISMATCH", cohort_value), call. = FALSE)
  }
  checks <- checks + 10L
}
cat("PATIENT_BALANCE_AGREEMENT_INDEPENDENT_QA_PASS checks=", checks, "\n", sep = "")
