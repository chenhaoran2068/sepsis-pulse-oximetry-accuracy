#!/usr/bin/env Rscript
# Independent numeric and file-contract check for the candidate BA renderer.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: MANIFEST.tsv WINDOW(60|5) FIGURE.pdf RECEIPT.tsv", call. = FALSE)
manifest <- read.delim(args[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
window <- suppressWarnings(as.integer(args[[2L]]))
figure <- args[[3L]]
receipt <- read.delim(args[[4L]], check.names = FALSE, stringsAsFactors = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (is.na(window) || !(window %in% c(5L, 60L))) stop("WINDOW_INVALID", call. = FALSE)
if (!identical(sort(manifest$cohort), sort(cohorts)) || !identical(sort(receipt$cohort), sort(cohorts))) {
  stop("COHORT_SET_INVALID", call. = FALSE)
}
if (!file.exists(figure) || file.info(figure)$size < 5000L) stop("FIGURE_EMPTY_OR_MISSING", call. = FALSE)
con <- file(figure, "rb")
magic <- readChar(con, 4L, useBytes = TRUE)
close(con)
if (magic != "%PDF") stop("FIGURE_NOT_PDF", call. = FALSE)

checks <- 0L
for (cohort in cohorts) {
  m <- manifest[manifest$cohort == cohort, , drop = FALSE]
  r <- receipt[receipt$cohort == cohort, , drop = FALSE]
  if (nrow(m) != 1L || nrow(r) != 1L) stop("COHORT_CARDINALITY_INVALID", call. = FALSE)
  x <- readRDS(m$analysis_rds[[1L]])
  s <- read.delim(m$agreement_tsv[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
  recomputed <- mean(as.numeric(x$spo2_saturation_percent) - as.numeric(x$sao2_saturation_percent))
  if (r$window_minutes[[1L]] != window || r$pair_n[[1L]] != nrow(x) ||
      abs(r$recomputed_mean_bias[[1L]] - recomputed) > 1e-10 ||
      abs(r$reported_mean_bias[[1L]] - s$mean_bias[[1L]]) > 1e-10 ||
      abs(r$lower_loa[[1L]] - s$lower_loa[[1L]]) > 1e-10 ||
      abs(r$upper_loa[[1L]] - s$upper_loa[[1L]]) > 1e-10) {
    stop(paste("RECEIPT_NOT_REPRODUCED", cohort), call. = FALSE)
  }
  checks <- checks + 6L
}
cat("BLAND_ALTMAN_INDEPENDENT_QA_PASS checks=", checks, " window=", window, "\n", sep = "")
