#!/usr/bin/env Rscript
# Independent per-term reporting-scale check; does not source plot renderer.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST.tsv WINDOW OUTPUT.pdf CHART_DATA.tsv RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
window <- as.integer(args[[2L]])
chart <- fread(args[[4L]])
receipt <- fread(args[[5L]])
checks <- character()
check <- function(name, ok) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(paste("QA_FAIL", name), call. = FALSE)
  checks <<- c(checks, name)
}
check("pdf_magic", identical(readBin(args[[3L]], "raw", n = 4L), charToRaw("%PDF")))
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
check("exact_grid", nrow(chart) == 50L && setequal(chart$cohort, cohorts) &&
      all(chart$window_minutes == window) && nrow(receipt) == 5L)
check("ten_rows_each", all(receipt$variable_n == 10L) &&
      all(receipt$estimated_n + receipt$categorical_age_n + receipt$not_in_model_n == 10L))
for (cohort_name in cohorts) {
  source <- fread(m[cohort == cohort_name, lmm_tsv])
  current <- chart[cohort == cohort_name]
  check(paste0(cohort_name, "_rows"), nrow(current) == 10L &&
        current[variable_id == "paired_mean", status] == "estimated")
  for (i in seq_len(nrow(current))) {
    z <- current[i]
    match_source <- source[term == z$term[[1L]]]
    if (z$status[[1L]] == "estimated") {
      check(paste0(cohort_name, "_term_unique_", i), nrow(match_source) == 1L)
      scale <- match_source$reporting_multiplier[[1L]]
      check(paste0(cohort_name, "_scaled_", i),
            abs(z$estimate_scaled[[1L]] - match_source$estimate[[1L]] * scale) < 1e-10 &&
            abs(z$lower_scaled[[1L]] - match_source$ci_lower[[1L]] * scale) < 1e-10 &&
            abs(z$upper_scaled[[1L]] - match_source$ci_upper[[1L]] * scale) < 1e-10)
      check(paste0(cohort_name, "_ci_", i),
            z$lower_scaled[[1L]] <= z$estimate_scaled[[1L]] &&
            z$estimate_scaled[[1L]] <= z$upper_scaled[[1L]])
    } else {
      check(paste0(cohort_name, "_omission_", i), nrow(match_source) == 0L &&
            is.na(z$estimate_scaled[[1L]]))
    }
  }
}
cat("LMM_FOREST_INDEPENDENT_QA_PASS checks=", length(checks), "\n", sep = "")
