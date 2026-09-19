#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: LMM60_MANIFEST.tsv LMM5_MANIFEST.tsv ROWS.tsv RECEIPT.tsv", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
rows <- fread(args[[3L]])
receipt <- fread(args[[4L]])
if (nrow(rows) != receipt$row_n[[1L]] || receipt$panel_n[[1L]] != 10L) {
  stop("LMM_BALANCE_OUTPUT_COUNT_INVALID", call. = FALSE)
}
checks <- 0L
for (k in seq_along(c(60L, 5L))) {
  window_value <- c(60L, 5L)[[k]]
  m <- fread(args[[k]])
  for (cohort_value in cohorts) {
    path <- m[cohort == cohort_value]$lmm_tsv
    p <- fread(path)
    b <- fread(file.path(dirname(path), "lmm_patient_balanced_pooled.tsv"))
    b <- b[match(p$term, b$term)]
    shown <- rows[cohort == unname(pretty[[cohort_value]]) & window_minutes == window_value]
    if (nrow(shown) != nrow(p)) stop("LMM_BALANCE_PANEL_SIZE_INVALID", call. = FALSE)
    for (i in seq_len(nrow(p))) {
      fmt <- function(z) sprintf("%+.2f (%+.2f to %+.2f)", z$estimate[[i]] * z$reporting_multiplier[[i]],
        z$ci_lower[[i]] * z$reporting_multiplier[[i]], z$ci_upper[[i]] * z$reporting_multiplier[[i]])
      if (shown$term[[i]] != p$term[[i]] || shown$label[[i]] != p$reporting_label[[i]] ||
          shown$primary_with_95ci[[i]] != fmt(p) || shown$balanced_with_95ci[[i]] != fmt(b) ||
          abs(shown$primary_estimate[[i]] - p$estimate[[i]] * p$reporting_multiplier[[i]]) > 1e-12 ||
          abs(shown$balanced_estimate[[i]] - b$estimate[[i]] * b$reporting_multiplier[[i]]) > 1e-12 ||
          shown$same_direction[[i]] != ifelse(sign(p$estimate[[i]]) == sign(b$estimate[[i]]), "Yes", "No")) {
        stop(paste("LMM_BALANCE_ROW_MISMATCH", cohort_value, window_value, p$term[[i]]), call. = FALSE)
      }
      checks <- checks + 7L
    }
  }
}
cat("LMM_BALANCE_TABLE_INDEPENDENT_QA_PASS checks=", checks, "\n", sep = "")
