#!/usr/bin/env Rscript
# Independent input-grid and receipt checks for the stratified agreement plot.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: MANIFEST.tsv OUTPUT.pdf RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
receipt <- fread(args[[3L]])
checks <- character()
check <- function(name, ok) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(paste("QA_FAIL", name), call. = FALSE)
  checks <<- c(checks, name)
}
check("pdf_magic", identical(readBin(args[[2L]], "raw", n = 4L), charToRaw("%PDF")))
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
check("receipt_cohort_set", nrow(receipt) == 5L && setequal(receipt$cohort, cohorts))
for (cohort_name in cohorts) {
  source <- file.path(dirname(m[cohort == cohort_name, agreement_tsv]),
                      "agreement_by_spo2_stratum.tsv")
  x <- fread(source)
  r <- receipt[cohort == cohort_name]
  check(paste0(cohort_name, "_grid"), nrow(x) == 4L &&
        setequal(x$spo2_stratum, strata) && all(x$window_minutes == 60L))
  check(paste0(cohort_name, "_receipt"), nrow(r) == 1L && r$stratum_n[[1L]] == 4L &&
        r$window_minutes[[1L]] == 60L && r$pair_n_total[[1L]] == sum(x$pair_n))
  for (i in seq_len(nrow(x))) {
    check(paste0(cohort_name, "_loa_", i),
          x$lower_loa[[i]] < x$mean_bias[[i]] && x$mean_bias[[i]] < x$upper_loa[[i]])
  }
}
cat("STRATUM_AGREEMENT_INDEPENDENT_QA_PASS checks=", length(checks), "\n", sep = "")
