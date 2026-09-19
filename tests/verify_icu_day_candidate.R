#!/usr/bin/env Rscript
# Independent aggregate/figure contract review for the candidate ICU-day plot.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: MANIFEST.tsv WINDOW(60|5) FIGURE.pdf RECEIPT.tsv", call. = FALSE)
manifest <- read.delim(args[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
window <- suppressWarnings(as.integer(args[[2L]]))
pdf_path <- args[[3L]]
receipt <- read.delim(args[[4L]], check.names = FALSE, stringsAsFactors = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
if (is.na(window) || !(window %in% c(5L, 60L))) stop("WINDOW_INVALID", call. = FALSE)
if (!setequal(manifest$cohort, cohorts) || !setequal(receipt$cohort, cohorts) ||
    nrow(manifest) != 5L || nrow(receipt) != 5L) stop("COHORT_SET_INVALID", call. = FALSE)
if (!file.exists(pdf_path) || file.info(pdf_path)$size < 5000L) stop("FIGURE_MISSING_OR_EMPTY", call. = FALSE)
con <- file(pdf_path, "rb")
magic <- readChar(con, 4L, useBytes = TRUE)
close(con)
if (magic != "%PDF") stop("FIGURE_NOT_PDF", call. = FALSE)
checks <- 0L
for (cohort in cohorts) {
  m <- manifest[manifest$cohort == cohort, , drop = FALSE]
  r <- receipt[receipt$cohort == cohort, , drop = FALSE]
  if (nrow(m) != 1L || nrow(r) != 1L) stop("COHORT_CARDINALITY_INVALID", call. = FALSE)
  x <- read.delim(m$day_tsv[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(x) != 28L || any(x$cohort != cohort) || any(x$window_minutes != window) ||
      !setequal(x$relative_icu_day, 1:7) || !setequal(x$spo2_stratum, strata) ||
      anyDuplicated(paste(x$relative_icu_day, x$spo2_stratum))) {
    stop(paste("SOURCE_GRID_INVALID", cohort), call. = FALSE)
  }
  if (r$cell_n[[1L]] != 28L ||
      r$displayable_n[[1L]] != sum(x$cell_status == "DISPLAYABLE") ||
      r$sparse_n[[1L]] != sum(x$cell_status != "DISPLAYABLE") ||
      r$window_minutes[[1L]] != window || r$lower_panel_max_pct[[1L]] < 10) {
    stop(paste("RECEIPT_NOT_REPRODUCED", cohort), call. = FALSE)
  }
  if (any(x$cell_status != "DISPLAYABLE" & x$denominator_n >= 20 & x$patient_n >= 10)) {
    stop(paste("SPARSE_RULE_INCONSISTENT", cohort), call. = FALSE)
  }
  checks <- checks + 7L
}
cat("ICU_DAY_INDEPENDENT_QA_PASS checks=", checks, " window=", window, "\n", sep = "")
