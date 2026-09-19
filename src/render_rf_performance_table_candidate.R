#!/usr/bin/env Rscript
# Technical candidate for Supplementary Table 23 from five cohort-specific forest outputs.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: RF_MANIFEST.tsv NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
if (!identical(names(m), c("cohort", "rf_dir")) || nrow(m) != 5L ||
    !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort)) {
  stop("RF_PERFORMANCE_MANIFEST_INVALID", call. = FALSE)
}
outputs <- args[2:4]
if (any(file.exists(outputs))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(outputs)))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
metrics <- c("rmse", "mae", "r_squared", "calibration_intercept", "calibration_slope")
rows <- rbindlist(lapply(cohorts, function(cohort_value) {
  dir <- m[cohort == cohort_value]$rf_dir
  perf_path <- file.path(dir, "rf_outer_test_performance.tsv")
  hyper_path <- file.path(dir, "rf_selected_hyperparameters.tsv")
  if (length(dir) != 1L || !file.exists(perf_path) || !file.exists(hyper_path)) {
    stop("RF_PERFORMANCE_SOURCE_MISSING", call. = FALSE)
  }
  p <- fread(perf_path)
  h <- fread(hyper_path)
  needed <- c("metric", "estimate", "ci_lower", "ci_upper", "bootstrap_requested",
              "bootstrap_successful", "cohort", "development_pair_n", "development_patient_n",
              "test_pair_n", "test_patient_n")
  if (!all(needed %in% names(p)) || nrow(p) != 5L ||
      !setequal(p$metric, metrics) || anyDuplicated(p$metric) ||
      !all(p$cohort == cohort_value) || nrow(h) != 1L || h$cohort[[1L]] != cohort_value ||
      h$resolved_mtry_min[[1L]] != h$resolved_mtry_max[[1L]] ||
      uniqueN(p$development_pair_n) != 1L || uniqueN(p$development_patient_n) != 1L ||
      uniqueN(p$test_pair_n) != 1L || uniqueN(p$test_patient_n) != 1L) {
    stop("RF_PERFORMANCE_SOURCE_CONTRACT_INVALID", call. = FALSE)
  }
  if (anyNA(p[, ..needed]) || any(!is.finite(as.matrix(p[, .(estimate, ci_lower, ci_upper)]))) ||
      any(p$ci_lower > p$estimate | p$estimate > p$ci_upper) ||
      any(p$bootstrap_requested < 1 | p$bootstrap_successful < 1 |
          p$bootstrap_successful > p$bootstrap_requested) ||
      any(p$development_patient_n < 1 | p$test_patient_n < 1 |
          p$development_pair_n < p$development_patient_n |
          p$test_pair_n < p$test_patient_n) ||
      h$resolved_mtry_min[[1L]] < 1 || h$min_node_size[[1L]] < 1) {
    stop("RF_PERFORMANCE_VALUES_INVALID", call. = FALSE)
  }
  p <- p[match(metrics, metric)]
  fmt <- function(row) sprintf("%.2f (%.2f to %.2f)", row$estimate, row$ci_lower, row$ci_upper)
  data.table(cohort = unname(pretty[[cohort_value]]), cohort_key = cohort_value,
    development_patients_pairs = paste0(p$development_patient_n[[1L]], "/", p$development_pair_n[[1L]]),
    heldout_patients_pairs = paste0(p$test_patient_n[[1L]], "/", p$test_pair_n[[1L]]),
    mtry_min_node_size = paste0(h$resolved_mtry_min[[1L]], "/", h$min_node_size[[1L]]),
    rmse_with_95ci = fmt(p[metric == "rmse"]),
    mae_with_95ci = fmt(p[metric == "mae"]),
    r_squared_with_95ci = fmt(p[metric == "r_squared"]),
    calibration_intercept_with_95ci = fmt(p[metric == "calibration_intercept"]),
    calibration_slope_with_95ci = fmt(p[metric == "calibration_slope"]),
    source_mode = trimws(readLines(file.path(dir, "status.txt"), warn = FALSE)[[1L]]))
}))
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(outputs[[1L]], width = 11.69, height = 8.27, family = "sans", bg = "white")
grid.newpage()
grid.text("Supplementary Table 23. Held-Out Random-Forest Performance",
          x = unit(0.42, "in"), y = unit(7.73, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 10.5, fontface = "bold"))
xs <- c(0.42, 1.85, 3.28, 4.55, 5.67, 6.82, 7.95, 9.03, 10.22)
headers <- c("Cohort", "Dev. patients/pairs", "Test patients/pairs", "mtry/node",
             "RMSE (95% CI)", "MAE (95% CI)", "R2 (95% CI)",
             "Cal. intercept (95% CI)", "Cal. slope (95% CI)")
columns <- c("cohort", "development_patients_pairs", "heldout_patients_pairs",
             "mtry_min_node_size", "rmse_with_95ci", "mae_with_95ci", "r_squared_with_95ci",
             "calibration_intercept_with_95ci", "calibration_slope_with_95ci")
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.85, 6.85), "in"), gp = gpar(lwd = 1.1))
for (j in seq_along(headers)) grid.text(headers[[j]], x = unit(xs[[j]], "in"),
  y = unit(6.58, "in"), just = c("left", "centre"), gp = gpar(fontsize = 6.7, fontface = "bold"))
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.32, 6.32), "in"), gp = gpar(lwd = 0.6))
for (i in seq_len(5L)) for (j in seq_along(columns)) {
  grid.text(as.character(rows[[columns[[j]]]][[i]]), x = unit(xs[[j]], "in"),
    y = unit(5.95 - (i - 1L) * 0.52, "in"), just = c("left", "centre"),
    gp = gpar(fontsize = 6.4))
}
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(3.54, 3.54), "in"), gp = gpar(lwd = 1.0))
grid.text("Technical candidate. Synthetic QA can use smoke and paper settings; mode is recorded in the rows file.",
  x = unit(0.42, "in"), y = unit(3.21, "in"), just = c("left", "centre"),
  gp = gpar(fontsize = 7.1))
dev.off()
fwrite(rows, outputs[[2L]], sep = "\t")
fwrite(data.table(table_id = "Supplementary Table 23", cohort_n = 5L, row_n = 5L,
  paper_setting_n = sum(rows$source_mode == "RF_PAPER_SETTINGS_COMPLETE_CANDIDATE_ONLY")),
  outputs[[3L]], sep = "\t")
cat("RF_PERFORMANCE_TABLE_CANDIDATE_PASS rows=5\n")
