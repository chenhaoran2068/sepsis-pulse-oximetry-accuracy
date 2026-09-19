#!/usr/bin/env Rscript
# Supplementary Tables 16 and 17 from independent threshold aggregates.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST.tsv KIND(reverse|patient) NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
kind <- args[[2L]]
pdf_path <- args[[3L]]
rows_path <- args[[4L]]
receipt_path <- args[[5L]]
if (!kind %in% c("reverse", "patient")) stop("KIND_INVALID", call. = FALSE)
if (!identical(names(m), c("cohort", "window_minutes", "low_sao2_tsv"))) stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
expected <- CJ(cohort = cohorts, window_minutes = c(5L, 60L))
actual <- unique(m[, .(cohort, window_minutes)])
if (nrow(m) != 10L || nrow(actual) != 10L ||
    nrow(fsetdiff(expected, actual)) || nrow(fsetdiff(actual, expected))) {
  stop("TEN_COHORT_WINDOW_SET_INVALID", call. = FALSE)
}
if (any(file.exists(c(pdf_path, rows_path, receipt_path)))) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(c(pdf_path, rows_path, receipt_path))))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
file_name <- if (kind == "reverse") "threshold_from_displayed_spo2.tsv" else "threshold_patient_summary.tsv"
names_pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
                  SICDB = "SICdb", Lianyungang = "Lianyungang")
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
source_rows <- lapply(seq_len(nrow(m)), function(i) {
  path <- file.path(dirname(m$low_sao2_tsv[[i]]), file_name)
  if (!file.exists(path)) stop(paste("AUX_THRESHOLD_INPUT_MISSING", m$cohort[[i]],
                                     m$window_minutes[[i]]), call. = FALSE)
  x <- fread(path)
  required <- if (kind == "reverse") {
    c("displayed_spo2_threshold_ge", "accepted_pair_n", "conditional_denominator_spo2_ge_n",
      "conditional_numerator_sao2_lt88_n", "conditional_proportion", "ci_lower", "ci_upper", "cell_status")
  } else {
    c("displayed_spo2_threshold_ge", "accepted_patient_n", "affected_patient_n",
      "affected_patient_proportion")
  }
  if (!all(required %in% names(x))) stop("AUX_THRESHOLD_COLUMNS_MISSING", call. = FALSE)
  if (nrow(x) != 2L || !setequal(x$displayed_spo2_threshold_ge, c(88, 92))) {
    stop("AUX_THRESHOLD_TWO_ROWS_INVALID", call. = FALSE)
  }
  x[, `:=`(cohort = m$cohort[[i]], window_minutes = m$window_minutes[[i]])]
  x
})
data <- rbindlist(source_rows, use.names = TRUE, fill = TRUE)
if (kind == "reverse") {
  if (anyNA(data[, .(accepted_pair_n, conditional_denominator_spo2_ge_n,
                      conditional_numerator_sao2_lt88_n, cell_status)]) ||
      any(!data$cell_status %in% c("ESTIMATED", "SPARSE")) ||
      any(data$conditional_denominator_spo2_ge_n < 1 |
          data$conditional_denominator_spo2_ge_n > data$accepted_pair_n |
          data$conditional_numerator_sao2_lt88_n < 0 |
          data$conditional_numerator_sao2_lt88_n > data$conditional_denominator_spo2_ge_n)) {
    stop("REVERSE_COUNTS_INVALID", call. = FALSE)
  }
  if (data[cell_status == "ESTIMATED", any(!is.finite(conditional_proportion) |
        !is.finite(ci_lower) | !is.finite(ci_upper) |
        abs(conditional_proportion -
        conditional_numerator_sao2_lt88_n / conditional_denominator_spo2_ge_n) > 1e-10 |
        ci_lower < 0 | ci_upper > 1 | ci_lower > conditional_proportion |
        ci_upper < conditional_proportion)]) stop("REVERSE_CI_INVALID", call. = FALSE)
  nesting <- data[, .(den88 = conditional_denominator_spo2_ge_n[displayed_spo2_threshold_ge == 88],
                      den92 = conditional_denominator_spo2_ge_n[displayed_spo2_threshold_ge == 92],
                      n88 = conditional_numerator_sao2_lt88_n[displayed_spo2_threshold_ge == 88],
                      n92 = conditional_numerator_sao2_lt88_n[displayed_spo2_threshold_ge == 92]),
                  by = .(cohort, window_minutes)]
  if (any(nesting$den92 > nesting$den88 | nesting$n92 > nesting$n88)) {
    stop("REVERSE_NESTING_INVALID", call. = FALSE)
  }
} else {
  if (anyNA(data[, .(accepted_patient_n, affected_patient_n, affected_patient_proportion)]) ||
      any(data$accepted_patient_n < 1 | data$affected_patient_n < 0 |
          data$affected_patient_n > data$accepted_patient_n |
          abs(data$affected_patient_proportion -
              data$affected_patient_n / data$accepted_patient_n) > 1e-10)) {
    stop("PATIENT_SUMMARY_VALUES_INVALID", call. = FALSE)
  }
  nesting <- data[, .(den_n = uniqueN(accepted_patient_n),
                      n88 = affected_patient_n[displayed_spo2_threshold_ge == 88],
                      n92 = affected_patient_n[displayed_spo2_threshold_ge == 92]),
                  by = .(cohort, window_minutes)]
  if (any(nesting$den_n != 1L | nesting$n92 > nesting$n88)) {
    stop("PATIENT_SUMMARY_NESTING_INVALID", call. = FALSE)
  }
}
rows <- rbindlist(lapply(seq_along(cohorts), function(i) {
  rbindlist(lapply(c(60L, 5L), function(w) {
    x <- data[cohort == cohorts[[i]] & window_minutes == w]
    a <- x[displayed_spo2_threshold_ge == 88]
    b <- x[displayed_spo2_threshold_ge == 92]
    if (kind == "reverse") {
      display <- function(z) {
        base <- paste0(fmt_n(z$conditional_numerator_sao2_lt88_n), "/",
                       fmt_n(z$conditional_denominator_spo2_ge_n))
        if (z$cell_status == "SPARSE") base else sprintf("%s (%.1f)", base,
                                                        100 * z$conditional_proportion)
      }
      ci <- function(z) if (z$cell_status == "SPARSE") "-" else
        sprintf("%.1f to %.1f", 100 * z$ci_lower, 100 * z$ci_upper)
      data.table(cohort = if (w == 60L) unname(names_pretty[[cohorts[[i]]]]) else "",
                 window = paste0(w, " min"), spo2_ge88_n_N_pct = display(a),
                 spo2_ge88_ci_pct = ci(a), spo2_ge92_n_N_pct = display(b),
                 spo2_ge92_ci_pct = ci(b))
    } else {
      display <- function(z) sprintf("%s/%s (%.1f)", fmt_n(z$affected_patient_n),
                                     fmt_n(z$accepted_patient_n), 100 * z$affected_patient_proportion)
      data.table(cohort = if (w == 60L) unname(names_pretty[[cohorts[[i]]]]) else "",
                 window = paste0(w, " min"), occult_patient_n_N_pct = display(a),
                 severe_patient_n_N_pct = display(b))
    }
  }))
}))
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
grid.newpage()
title <- if (kind == "reverse") {
  "Supplementary Table 16. Arterial Hypoxemia Among Pairs Above Displayed SpO2 Thresholds"
} else {
  "Supplementary Table 17. Patients With an Observed Occult Hypoxemia Episode"
}
grid.text(title, x = unit(0.42, "in"), y = unit(7.72, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 11, fontface = "bold"))
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.88, 6.88), "in"), gp = gpar(lwd = 1.2))
if (kind == "reverse") {
  xs <- c(0.48, 1.85, 2.88, 5.28, 6.78, 9.20)
  headers <- c("Cohort", "Window", "SpO2 >=88%: SaO2 <88%, n/N (%)", "95% CI, %",
               "SpO2 >=92%: SaO2 <88%, n/N (%)", "95% CI, %")
} else {
  xs <- c(0.48, 2.55, 4.12, 7.66)
  headers <- c("Cohort", "Window", "At least one occult episode, n/N (%)",
               "At least one severe occult episode, n/N (%)")
}
for (j in seq_along(headers)) {
  grid.text(headers[[j]], x = unit(xs[[j]], "in"), y = unit(6.56, "in"),
            just = c("left", "centre"), gp = gpar(fontsize = 8.0, fontface = "bold"))
}
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.29, 6.29), "in"), gp = gpar(lwd = 0.7))
ys <- seq(6.02, 2.42, by = -0.40)
for (i in seq_len(nrow(rows))) {
  if (i > 1L && i %% 2L == 1L) {
    yy <- (ys[[i - 1L]] + ys[[i]]) / 2
    grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(yy, yy), "in"),
               gp = gpar(lwd = 0.35, col = "grey65"))
  }
  for (j in seq_len(ncol(rows))) {
    grid.text(as.character(rows[[j]][[i]]), x = unit(xs[[j]], "in"), y = unit(ys[[i]], "in"),
              just = c("left", "centre"), gp = gpar(fontsize = 8.0))
  }
}
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(2.18, 2.18), "in"), gp = gpar(lwd = 1.1))
notes <- if (kind == "reverse") {
  c("The denominator is eligible pairs at or above each displayed SpO2 threshold. The numerator has SaO2 <88%.",
    "Confidence intervals use patient-clustered bootstrap resamples.",
    "Abbreviations: CI, confidence interval. SaO2, arterial oxygen saturation. SpO2, peripheral oxygen saturation.")
} else {
  c("Each patient is counted once for each episode pattern, regardless of their number of qualifying pairs.",
    "The denominator is all patients contributing an eligible pair to the corresponding pairing analysis.",
    "Abbreviations: SaO2, arterial oxygen saturation. SpO2, peripheral oxygen saturation.")
}
for (i in seq_along(notes)) {
  grid.text(notes[[i]], x = unit(0.43, "in"), y = unit(1.88 - (i - 1L) * 0.29, "in"),
            just = c("left", "centre"), gp = gpar(fontsize = 7.2))
}
dev.off()
fwrite(rows, rows_path, sep = "\t")
fwrite(data.table(table_id = if (kind == "reverse") "Supplementary Table 16" else "Supplementary Table 17",
                  cohort_n = 5L, analysis_n = 10L, input_row_n = nrow(data), display_row_n = nrow(rows)),
       receipt_path, sep = "\t")
cat("THRESHOLD_AUX_TABLE_CANDIDATE_RENDER_PASS kind=", kind, " rows=10\n", sep = "")
