#!/usr/bin/env Rscript
# Technical candidate for Supplementary Table 13 from cohort/window agreement aggregates.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST60.tsv MANIFEST5.tsv NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
labels <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
pdf_path <- args[[3L]]
rows_path <- args[[4L]]
receipt_path <- args[[5L]]
outputs <- c(pdf_path, rows_path, receipt_path)
if (any(file.exists(outputs))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(outputs)))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
records <- lapply(seq_along(c(60L, 5L)), function(k) {
  window <- c(60L, 5L)[[k]]
  m <- fread(args[[k]])
  if (!identical(names(m), c("cohort", "analysis_rds", "agreement_tsv")) ||
      nrow(m) != 5L || !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort)) {
    stop("STRATUM_MANIFEST_INVALID", call. = FALSE)
  }
  lapply(seq_len(nrow(m)), function(i) {
    overall_path <- m$agreement_tsv[[i]]
    strata_path <- file.path(dirname(overall_path), "agreement_by_spo2_stratum.tsv")
    if (!file.exists(overall_path) || !file.exists(strata_path)) {
      stop("STRATUM_AGREEMENT_INPUT_MISSING", call. = FALSE)
    }
    overall <- fread(overall_path)
    z <- fread(strata_path)
    needed <- c("cohort", "window_minutes", "spo2_stratum", "pair_n", "patient_n",
                "mean_bias", "mean_bias_ci_lower", "mean_bias_ci_upper",
                "lower_loa", "lower_loa_ci_lower", "lower_loa_ci_upper",
                "upper_loa", "upper_loa_ci_lower", "upper_loa_ci_upper")
    if (nrow(overall) != 1L || nrow(z) != 4L || !all(needed %in% names(z)) ||
        !setequal(z$spo2_stratum, strata) || anyDuplicated(z$spo2_stratum) ||
        !all(z$cohort == m$cohort[[i]]) || !all(z$window_minutes == window) ||
        overall$cohort[[1L]] != m$cohort[[i]] || overall$window_minutes[[1L]] != window ||
        sum(z$pair_n) != overall$pair_n[[1L]]) {
      stop("STRATUM_AGREEMENT_CONTRACT_INVALID", call. = FALSE)
    }
    numeric_fields <- setdiff(needed, c("cohort", "spo2_stratum"))
    if (anyNA(z[, ..numeric_fields]) || any(!is.finite(as.matrix(z[, ..numeric_fields]))) ||
        any(z$pair_n < 1 | z$patient_n < 1 | z$patient_n > z$pair_n) ||
        any(z$mean_bias_ci_lower > z$mean_bias | z$mean_bias > z$mean_bias_ci_upper |
            z$lower_loa_ci_lower > z$lower_loa | z$lower_loa > z$lower_loa_ci_upper |
            z$upper_loa_ci_lower > z$upper_loa | z$upper_loa > z$upper_loa_ci_upper |
            z$lower_loa > z$mean_bias | z$mean_bias > z$upper_loa)) {
      stop("STRATUM_AGREEMENT_VALUES_INVALID", call. = FALSE)
    }
    z
  })
})
data <- rbindlist(unlist(records, recursive = FALSE), use.names = TRUE, fill = TRUE)
fmt_ci <- function(x, lo, hi) sprintf("%.2f (%.2f to %.2f)", x, lo, hi)
rows <- rbindlist(lapply(c(60L, 5L), function(window_value) {
  rbindlist(lapply(cohorts, function(cohort_value) {
    rbindlist(lapply(seq_along(strata), function(j) {
      z <- data[cohort == cohort_value & window_minutes == window_value & spo2_stratum == strata[[j]]]
      if (nrow(z) != 1L) stop("STRATUM_LOOKUP_NOT_UNIQUE", call. = FALSE)
      data.table(cohort = unname(pretty[[cohort_value]]), window = paste0(window_value, " min"),
                 spo2_stratum = labels[[j]], pairs_n = format(z$pair_n, big.mark = ",", trim = TRUE),
                 mean_bias_with_95ci = fmt_ci(z$mean_bias, z$mean_bias_ci_lower, z$mean_bias_ci_upper),
                 lower_loa_with_95ci = fmt_ci(z$lower_loa, z$lower_loa_ci_lower, z$lower_loa_ci_upper),
                 upper_loa_with_95ci = fmt_ci(z$upper_loa, z$upper_loa_ci_lower, z$upper_loa_ci_upper))
    }))
  }))
}))
if (nrow(rows) != 40L) stop("STRATUM_OUTPUT_ROW_COUNT_INVALID", call. = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
xs <- c(0.42, 2.22, 3.38, 4.55, 6.7, 8.9)
headers <- c("Cohort", "SpO2 stratum", "Pairs, n", "Mean Bias (95% CI)",
             "Lower LoA (95% CI)", "Upper LoA (95% CI)")
for (window_value in c(60L, 5L)) {
  grid.newpage()
  grid.text("Supplementary Table 13. SpO2 Stratum-Specific Agreement Results",
            x = unit(0.42, "in"), y = unit(7.75, "in"), just = c("left", "centre"),
            gp = gpar(fontsize = 10.5, fontface = "bold"))
  grid.text(if (window_value == 60L) "60-minute primary analysis" else "5-minute sensitivity analysis",
            x = unit(0.42, "in"), y = unit(7.25, "in"), just = c("left", "centre"),
            gp = gpar(fontsize = 9, fontface = "bold"))
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.98, 6.98), "in"), gp = gpar(lwd = 1.1))
  for (j in seq_along(headers)) grid.text(headers[[j]], x = unit(xs[[j]], "in"),
      y = unit(6.75, "in"), just = c("left", "centre"), gp = gpar(fontsize = 7.8, fontface = "bold"))
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.51, 6.51), "in"), gp = gpar(lwd = 0.6))
  page <- rows[window == paste0(window_value, " min")]
  ys <- seq(6.28, by = -0.255, length.out = 20L)
  for (i in seq_len(nrow(page))) {
    if (i %in% c(5L, 9L, 13L, 17L)) grid.lines(x = unit(c(0.42, 11.27), "in"),
      y = unit(c(ys[[i]] + 0.13, ys[[i]] + 0.13), "in"), gp = gpar(lwd = 0.35, col = "grey60"))
    for (j in seq_len(ncol(page))) {
      if (j == 2L) next
      column <- if (j > 2L) j - 1L else j
      grid.text(as.character(page[[j]][[i]]), x = unit(xs[[column]], "in"),
                y = unit(ys[[i]], "in"), just = c("left", "centre"), gp = gpar(fontsize = 7.1))
    }
  }
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(1.31, 1.31), "in"), gp = gpar(lwd = 1.0))
  grid.text("Bias = SpO2 - SaO2. LoA = 95% limits of agreement. CI = confidence interval.",
            x = unit(0.42, "in"), y = unit(1.04, "in"), just = c("left", "centre"),
            gp = gpar(fontsize = 7.2))
}
dev.off()
fwrite(rows, rows_path, sep = "\t")
fwrite(data.table(table_id = "Supplementary Table 13", cohort_n = 5L, analysis_n = 10L,
                  stratum_n = 4L, row_n = nrow(rows)), receipt_path, sep = "\t")
cat("STRATIFIED_AGREEMENT_TABLE_CANDIDATE_PASS rows=40\n")
