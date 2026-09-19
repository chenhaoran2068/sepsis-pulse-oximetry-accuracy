#!/usr/bin/env Rscript
# Read display cells independently from producer aggregates; no renderer sourcing.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST.tsv WINDOW OUTPUT.pdf ROWS.tsv RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
window <- as.integer(args[[2L]])
pdf_path <- args[[3L]]
actual <- fread(args[[4L]])
receipt <- fread(args[[5L]])
checks <- character()
check <- function(name, ok) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(paste("QA_FAIL", name), call. = FALSE)
  checks <<- c(checks, name)
}
check("pdf_magic", identical(readBin(pdf_path, "raw", n = 4L), charToRaw("%PDF")))
check("five_rows", nrow(actual) == 5L && ncol(actual) == 6L)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
labels <- c("MIMIC", "AmsterdamUMCdb", "eICU", "SICdb", "Lianyungang")
check("order", identical(actual$cohort, labels))
fmt_ci <- function(z, stem) sprintf("%.2f (%.2f to %.2f)", z[[stem]],
                                  z[[paste0(stem, "_ci_lower")]],
                                  z[[paste0(stem, "_ci_upper")]])
pair_total <- 0L
for (i in seq_along(cohorts)) {
  source <- fread(m[cohort == cohorts[[i]], agreement_tsv])
  check(paste0(cohorts[[i]], "_source_row"), nrow(source) == 1L && source$window_minutes[[1L]] == window)
  pair_total <- pair_total + source$pair_n[[1L]]
  expected_count <- sprintf("%s (%s)", format(source$pair_n, big.mark = ",", scientific = FALSE,
                                              trim = TRUE),
                            format(source$stay_n, big.mark = ",", scientific = FALSE, trim = TRUE))
  check(paste0(cohorts[[i]], "_counts"), identical(actual$pairs_and_stays[[i]], expected_count))
  for (stem in c("mean_bias", "lower_loa", "upper_loa", "arms")) {
    check(paste0(cohorts[[i]], "_", stem),
          identical(actual[[paste0(stem, "_ci")]][[i]], fmt_ci(source, stem)))
  }
}
check("receipt_window", nrow(receipt) == 1L && receipt$window_minutes[[1L]] == window)
check("receipt_counts", receipt$row_n[[1L]] == 5L && receipt$cohort_n[[1L]] == 5L &&
      receipt$pair_n_total[[1L]] == pair_total)
cat("AGREEMENT_TABLE_INDEPENDENT_QA_PASS checks=", length(checks), "\n", sep = "")
