#!/usr/bin/env Rscript
# Independently compare every candidate cell with the source aggregates.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST.tsv KIND OUTPUT.pdf ROWS.tsv RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
kind <- args[[2L]]
rows <- fread(args[[4L]], colClasses = "character")
receipt <- fread(args[[5L]])
checks <- character()
check <- function(name, ok) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(paste("QA_FAIL", name), call. = FALSE)
  checks <<- c(checks, name)
}
check("pdf_magic", identical(readBin(args[[3L]], "raw", n = 4L), charToRaw("%PDF")))
check("shape", nrow(rows) == 10L && ncol(rows) == if (kind == "reverse") 6L else 4L)
check("receipt", nrow(receipt) == 1L && receipt$cohort_n[[1L]] == 5L &&
      receipt$analysis_n[[1L]] == 10L && receipt$input_row_n[[1L]] == 20L)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
labels <- c("MIMIC", "AmsterdamUMCdb", "eICU", "SICdb", "Lianyungang")
fmt_n <- function(z) format(z, big.mark = ",", scientific = FALSE, trim = TRUE)
for (i in seq_along(cohorts)) for (window in c(60L, 5L)) {
  parent <- dirname(m[cohort == cohorts[[i]] & window_minutes == window, low_sao2_tsv])
  source <- fread(file.path(parent, if (kind == "reverse")
    "threshold_from_displayed_spo2.tsv" else "threshold_patient_summary.tsv"))
  row_index <- 2L * (i - 1L) + if (window == 60L) 1L else 2L
  row <- rows[row_index]
  check(paste0("row_label_", i, "_", window),
        identical(row$cohort[[1L]], if (window == 60L) labels[[i]] else "") &&
        identical(row$window[[1L]], paste0(window, " min")))
  for (cutoff in c(88L, 92L)) {
    x <- source[displayed_spo2_threshold_ge == cutoff]
    check(paste0("input_cutoff_", i, "_", window, "_", cutoff), nrow(x) == 1L)
    field <- if (cutoff == 88L) "spo2_ge88" else "spo2_ge92"
    if (kind == "reverse") {
      expected_n <- sprintf("%s/%s (%.1f)", fmt_n(x$conditional_numerator_sao2_lt88_n),
                            fmt_n(x$conditional_denominator_spo2_ge_n), 100 * x$conditional_proportion)
      expected_ci <- sprintf("%.1f to %.1f", 100 * x$ci_lower, 100 * x$ci_upper)
      check(paste0("cell_n_", i, "_", window, "_", cutoff),
            identical(row[[paste0(field, "_n_N_pct")]][[1L]], expected_n))
      check(paste0("cell_ci_", i, "_", window, "_", cutoff),
            identical(row[[paste0(field, "_ci_pct")]][[1L]], expected_ci))
    } else {
      expected <- sprintf("%s/%s (%.1f)", fmt_n(x$affected_patient_n),
                          fmt_n(x$accepted_patient_n), 100 * x$affected_patient_proportion)
      field <- if (cutoff == 88L) "occult_patient_n_N_pct" else "severe_patient_n_N_pct"
      check(paste0("cell_patient_", i, "_", window, "_", cutoff),
            identical(row[[field]][[1L]], expected))
    }
  }
}
cat("THRESHOLD_AUX_INDEPENDENT_QA_PASS kind=", kind,
    " checks=", length(checks), "\n", sep = "")
