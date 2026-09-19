#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: MANIFEST60.tsv RECEIPT.tsv", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
m <- fread(args[[1L]])
r <- fread(args[[2L]])
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
if (nrow(r) != 20L) stop("STRATIFIED_BA_RECEIPT_ROW_COUNT_INVALID", call. = FALSE)
checks <- 0L
for (cohort_value in cohorts) {
  row <- m[cohort == cohort_value]
  if (nrow(row) != 1L) stop("SOURCE_MANIFEST_COHORT_MISSING", call. = FALSE)
  x <- as.data.table(readRDS(row$analysis_rds))
  x[, bias := spo2_saturation_percent - sao2_saturation_percent]
  x[, stratum := fifelse(spo2_saturation_percent < 88, strata[[1L]],
    fifelse(spo2_saturation_percent < 92, strata[[2L]],
      fifelse(spo2_saturation_percent <= 96, strata[[3L]], strata[[4L]])))]
  source <- fread(file.path(dirname(row$agreement_tsv), "agreement_by_spo2_stratum.tsv"))
  for (stratum_value in strata) {
    shown <- r[cohort == cohort_value & spo2_stratum == stratum_value]
    original <- source[cohort == cohort_value & spo2_stratum == stratum_value]
    subset <- x[stratum == stratum_value]
    if (nrow(shown) != 1L || nrow(original) != 1L ||
        shown$pair_n != nrow(subset) || shown$pair_n != original$pair_n ||
        abs(shown$mean_bias - original$mean_bias) > 1e-12 ||
        abs(shown$lower_loa - original$lower_loa) > 1e-12 ||
        abs(shown$upper_loa - original$upper_loa) > 1e-12) {
      stop(paste("STRATIFIED_BA_RECEIPT_MISMATCH", cohort_value, stratum_value), call. = FALSE)
    }
    checks <- checks + 6L
  }
}
cat("STRATIFIED_BA_INDEPENDENT_QA_PASS checks=", checks, "\n", sep = "")
