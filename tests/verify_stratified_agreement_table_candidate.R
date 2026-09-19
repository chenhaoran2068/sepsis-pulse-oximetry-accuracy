#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: MANIFEST60.tsv MANIFEST5.tsv ROWS.tsv RECEIPT.tsv", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
pretty <- c("MIMIC", "AmsterdamUMCdb", "eICU", "SICdb", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
actual <- fread(args[[3L]])
receipt <- fread(args[[4L]])
if (nrow(actual) != 40L || receipt$row_n[[1L]] != 40L) stop("STRATUM_OUTPUT_COUNT_INVALID", call. = FALSE)
format_three <- function(z, prefix) {
  sprintf("%.2f (%.2f to %.2f)", z[[prefix]], z[[paste0(prefix, "_ci_lower")]],
          z[[paste0(prefix, "_ci_upper")]])
}
checks <- 0L
for (window_index in seq_along(c(60L, 5L))) {
  window <- c(60L, 5L)[[window_index]]
  manifest <- fread(args[[window_index]])
  for (cohort_index in seq_along(cohorts)) {
    row <- manifest[cohort == cohorts[[cohort_index]]]
    if (nrow(row) != 1L) stop("SOURCE_MANIFEST_COHORT_MISSING", call. = FALSE)
    source <- fread(file.path(dirname(row$agreement_tsv), "agreement_by_spo2_stratum.tsv"))
    overall <- fread(row$agreement_tsv)
    if (sum(source$pair_n) != overall$pair_n[[1L]]) stop("SOURCE_STRATA_DO_NOT_SUM_TO_OVERALL", call. = FALSE)
    for (stratum_index in seq_along(strata)) {
      z <- source[spo2_stratum == strata[[stratum_index]]]
      position <- (window_index - 1L) * 20L + (cohort_index - 1L) * 4L + stratum_index
      shown <- actual[position]
      expected_pairs <- format(z$pair_n, big.mark = ",", trim = TRUE)
      if (shown$cohort != pretty[[cohort_index]] || shown$window != paste0(window, " min") ||
          shown$spo2_stratum != strata[[stratum_index]] || shown$pairs_n != expected_pairs ||
          shown$mean_bias_with_95ci != format_three(z, "mean_bias") ||
          shown$lower_loa_with_95ci != format_three(z, "lower_loa") ||
          shown$upper_loa_with_95ci != format_three(z, "upper_loa")) {
        stop(paste("STRATUM_DISPLAY_MISMATCH", window, cohorts[[cohort_index]], strata[[stratum_index]]), call. = FALSE)
      }
      checks <- checks + 7L
    }
  }
}
cat("STRATIFIED_AGREEMENT_TABLE_INDEPENDENT_QA_PASS checks=", checks, "\n", sep = "")
