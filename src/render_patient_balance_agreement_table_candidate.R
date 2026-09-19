#!/usr/bin/env Rscript
# Technical candidate for Supplementary Table 12 from 60-minute agreement outputs.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: MANIFEST.tsv NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
manifest <- fread(args[[1L]])
pdf_path <- args[[2L]]
rows_path <- args[[3L]]
receipt_path <- args[[4L]]
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
if (!identical(names(manifest), c("cohort", "analysis_rds", "agreement_tsv")) ||
    nrow(manifest) != 5L || !setequal(manifest$cohort, cohorts) || anyDuplicated(manifest$cohort)) {
  stop("BALANCE_AGREEMENT_MANIFEST_INVALID", call. = FALSE)
}
outputs <- c(pdf_path, rows_path, receipt_path)
if (any(file.exists(outputs))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(outputs)))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
source_rows <- lapply(seq_len(nrow(manifest)), function(i) {
  overall_path <- manifest$agreement_tsv[[i]]
  balance_path <- file.path(dirname(overall_path), "agreement_patient_balanced.tsv")
  if (!file.exists(overall_path) || !file.exists(balance_path)) stop("BALANCE_AGREEMENT_INPUT_MISSING", call. = FALSE)
  o <- fread(overall_path)
  b <- fread(balance_path)
  o_fields <- c("cohort", "window_minutes", "patient_n", "pair_n", "between_patient_variance",
                "between_stay_within_patient_variance", "within_stay_residual_variance",
                "total_difference_sd", "empirical_p2_5", "empirical_p2_5_ci_lower",
                "empirical_p2_5_ci_upper", "empirical_p97_5", "empirical_p97_5_ci_lower",
                "empirical_p97_5_ci_upper")
  b_fields <- c("cohort", "window_minutes", "patient_n", "pair_n", "balanced_mean_bias",
                "balanced_mean_bias_ci_lower", "balanced_mean_bias_ci_upper",
                "balanced_lower_loa", "balanced_lower_loa_ci_lower", "balanced_lower_loa_ci_upper",
                "balanced_upper_loa", "balanced_upper_loa_ci_lower", "balanced_upper_loa_ci_upper",
                "balanced_arms", "balanced_arms_ci_lower", "balanced_arms_ci_upper")
  if (nrow(o) != 1L || nrow(b) != 1L || !all(o_fields %in% names(o)) ||
      !all(b_fields %in% names(b)) || o$cohort[[1L]] != manifest$cohort[[i]] ||
      b$cohort[[1L]] != manifest$cohort[[i]] || o$window_minutes[[1L]] != 60L ||
      b$window_minutes[[1L]] != 60L || o$patient_n[[1L]] != b$patient_n[[1L]] ||
      o$pair_n[[1L]] != b$pair_n[[1L]]) {
    stop("BALANCE_AGREEMENT_SOURCE_CONTRACT_INVALID", call. = FALSE)
  }
  complete_o <- setdiff(o_fields, "between_stay_within_patient_variance")
  if (anyNA(o[, ..complete_o]) || anyNA(b[, ..b_fields]) ||
      any(!is.finite(as.numeric(unlist(o[, ..complete_o][, -c("cohort"), with = FALSE])))) ||
      any(!is.finite(as.numeric(unlist(b[, ..b_fields][, -c("cohort"), with = FALSE]))))) {
    stop("BALANCE_AGREEMENT_MISSING_OR_NONFINITE", call. = FALSE)
  }
  stay_variance <- o$between_stay_within_patient_variance[[1L]]
  if ((manifest$cohort[[i]] == "eICU") != is.finite(stay_variance)) {
    stop("BALANCE_AGREEMENT_STAY_COMPONENT_INVALID", call. = FALSE)
  }
  variance_sum <- o$between_patient_variance[[1L]] +
    (if (is.finite(stay_variance)) stay_variance else 0) +
    o$within_stay_residual_variance[[1L]]
  if (abs(sqrt(variance_sum) - o$total_difference_sd[[1L]]) > 1e-7) {
    stop("BALANCE_AGREEMENT_VARIANCE_IDENTITY_INVALID", call. = FALSE)
  }
  for (name in c("balanced_mean_bias", "balanced_lower_loa", "balanced_upper_loa", "balanced_arms")) {
    if (b[[paste0(name, "_ci_lower")]][[1L]] > b[[name]][[1L]] ||
        b[[name]][[1L]] > b[[paste0(name, "_ci_upper")]][[1L]]) {
      stop("BALANCE_AGREEMENT_CI_INVALID", call. = FALSE)
    }
  }
  list(overall = o, balanced = b)
})
names(source_rows) <- manifest$cohort
fmt <- function(value, lower, upper) sprintf("%.2f (%.2f to %.2f)", value, lower, upper)
fmt_column <- function(z, name) fmt(z[[name]], z[[paste0(name, "_ci_lower")]],
                                    z[[paste0(name, "_ci_upper")]])
rows <- rbindlist(lapply(cohorts, function(cohort_value) {
  o <- source_rows[[cohort_value]]$overall
  b <- source_rows[[cohort_value]]$balanced
  data.table(cohort = unname(pretty[[cohort_value]]),
    patients_n = format(o$patient_n, big.mark = ",", trim = TRUE),
    balanced_mean_bias_with_95ci = fmt_column(b, "balanced_mean_bias"),
    balanced_lower_loa_with_95ci = fmt_column(b, "balanced_lower_loa"),
    balanced_upper_loa_with_95ci = fmt_column(b, "balanced_upper_loa"),
    balanced_arms_with_95ci = fmt_column(b, "balanced_arms"),
    empirical_p2_5_with_95ci = fmt_column(o, "empirical_p2_5"),
    empirical_p97_5_with_95ci = fmt_column(o, "empirical_p97_5"),
    patient_variance = sprintf("%.2f", o$between_patient_variance),
    stay_within_patient_variance = if (is.finite(o$between_stay_within_patient_variance))
      sprintf("%.2f", o$between_stay_within_patient_variance) else "-",
    residual_variance = sprintf("%.2f", o$within_stay_residual_variance),
    total_difference_sd = sprintf("%.2f", o$total_difference_sd))
}))
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
grid.newpage()
grid.text("Supplementary Table 12. Patient-Balanced Agreement and Variance Components",
          x = unit(0.42, "in"), y = unit(7.73, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 10.5, fontface = "bold"))
draw_section <- function(heading, headers, columns, xs, top, row_ys) {
  grid.text(heading, x = unit(0.43, "in"), y = unit(top, "in"), just = c("left", "centre"),
            gp = gpar(fontsize = 9.0, fontface = "bold"))
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(top - 0.23, top - 0.23), "in"),
             gp = gpar(lwd = 1.05))
  for (j in seq_along(headers)) {
    grid.text(headers[[j]], x = unit(xs[[j]], "in"), y = unit(top - 0.48, "in"),
              just = c("left", "centre"), gp = gpar(fontsize = 7.1, fontface = "bold"))
  }
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(top - 0.72, top - 0.72), "in"),
             gp = gpar(lwd = 0.55))
  for (i in seq_len(5L)) for (j in seq_along(columns)) {
    grid.text(as.character(rows[[columns[[j]]]][[i]]), x = unit(xs[[j]], "in"),
              y = unit(row_ys[[i]], "in"), just = c("left", "centre"), gp = gpar(fontsize = 7.0))
  }
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(min(row_ys) - 0.27, min(row_ys) - 0.27), "in"),
             gp = gpar(lwd = 0.9))
}
draw_section("Patient-balanced sensitivity analysis",
  c("Cohort", "Patients", "Mean Bias (95% CI)", "Lower LoA (95% CI)",
    "Upper LoA (95% CI)", "A_rms (95% CI)"),
  c("cohort", "patients_n", "balanced_mean_bias_with_95ci", "balanced_lower_loa_with_95ci",
    "balanced_upper_loa_with_95ci", "balanced_arms_with_95ci"),
  c(0.43, 2.1, 3.05, 5.12, 7.18, 9.24), 7.20, seq(6.18, by = -0.37, length.out = 5L))
draw_section("Empirical distributional checks and variance components",
  c("Cohort", "P2.5 (95% CI)", "P97.5 (95% CI)", "Patient var.",
    "Stay var.", "Residual var.", "Total SD"),
  c("cohort", "empirical_p2_5_with_95ci", "empirical_p97_5_with_95ci",
    "patient_variance", "stay_within_patient_variance", "residual_variance", "total_difference_sd"),
  c(0.43, 2.1, 4.15, 6.3, 7.5, 8.7, 10.1), 4.05, seq(3.03, by = -0.37, length.out = 5L))
grid.text("Patient-balanced estimates give each patient total weight 1. Empirical quantiles are unweighted.",
          x = unit(0.43, "in"), y = unit(0.83, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 7.2))
dev.off()
fwrite(rows, rows_path, sep = "\t")
fwrite(data.table(table_id = "Supplementary Table 12", cohort_n = 5L, analysis_n = 5L,
                  row_n = nrow(rows)), receipt_path, sep = "\t")
cat("PATIENT_BALANCE_AGREEMENT_TABLE_CANDIDATE_PASS rows=5\n")
