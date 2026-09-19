#!/usr/bin/env Rscript
# Independent comparison of all rendered cells with source threshold aggregates.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: MANIFEST.tsv OUTPUT.pdf ROWS.tsv RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
actual <- fread(args[[3L]], colClasses = "character")
receipt <- fread(args[[4L]])
checks <- character()
check <- function(name, ok) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(paste("QA_FAIL", name), call. = FALSE)
  checks <<- c(checks, name)
}
check("pdf_magic", identical(readBin(args[[2L]], "raw", n = 4L), charToRaw("%PDF")))
check("row_shape", nrow(actual) == 10L && ncol(actual) == 7L)
check("receipt", nrow(receipt) == 1L && receipt$cohort_n[[1L]] == 5L &&
      receipt$analysis_n[[1L]] == 10L && receipt$input_row_n[[1L]] == 20L &&
      receipt$display_row_n[[1L]] == 10L)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
names_expected <- c("MIMIC", "AmsterdamUMCdb", "eICU", "SICdb", "Lianyungang")
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
for (i in seq_along(cohorts)) for (window in c(60L, 5L)) {
  source_path <- m[cohort == cohorts[[i]] & window_minutes == window, low_sao2_tsv]
  check(paste0("manifest_", i, "_", window), length(source_path) == 1L)
  source <- fread(source_path)
  check(paste0("input_", i, "_", window), nrow(source) == 2L &&
        setequal(source$threshold_spo2_ge, c(88L, 92L)))
  row_index <- 2L * (i - 1L) + if (window == 60L) 1L else 2L
  row <- actual[row_index]
  check(paste0("label_", i, "_", window),
        identical(row$cohort[[1L]], if (window == 60L) names_expected[[i]] else "") &&
        identical(row$window[[1L]], paste0(window, " min")))
  check(paste0("denom_", i, "_", window),
        identical(row$sao2_lt88_n[[1L]], fmt_n(source$pair_denominator_sao2_lt88_n[[1L]])))
  for (cutoff in c(88L, 92L)) {
    x <- source[threshold_spo2_ge == cutoff]
    prefix <- if (cutoff == 88L) "occult" else "severe"
    expected_episode <- sprintf("%s/%s (%.1f)", fmt_n(x$pair_numerator_n),
                                fmt_n(x$pair_denominator_sao2_lt88_n),
                                100 * x$pair_conditional_proportion)
    expected_ci <- sprintf("%.1f to %.1f", 100 * x$pair_ci_lower, 100 * x$pair_ci_upper)
    check(paste0(prefix, "_count_", i, "_", window),
          identical(row[[paste0(prefix, "_n_N_pct")]][[1L]], expected_episode))
    check(paste0(prefix, "_ci_", i, "_", window),
          identical(row[[paste0(prefix, "_ci_pct")]][[1L]], expected_ci))
  }
}
cat("THRESHOLD_TABLE_INDEPENDENT_QA_PASS checks=", length(checks), "\n", sep = "")
