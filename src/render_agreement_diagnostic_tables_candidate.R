#!/usr/bin/env Rscript
# Technical candidate for Supplementary Tables 14 and 15. Invented inputs only in tests.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST.tsv KIND(proportional|heteroscedastic) NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
manifest <- fread(args[[1L]])
kind <- args[[2L]]
pdf_path <- args[[3L]]
rows_path <- args[[4L]]
receipt_path <- args[[5L]]
if (!kind %in% c("proportional", "heteroscedastic")) stop("DIAGNOSTIC_KIND_INVALID", call. = FALSE)
if (!identical(names(manifest), c("cohort", "window_minutes", "diagnostic_dir"))) stop("DIAGNOSTIC_MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
expected <- CJ(cohort = cohorts, window_minutes = c(5L, 60L))
actual <- unique(manifest[, .(cohort, window_minutes)])
if (nrow(manifest) != 10L || nrow(actual) != 10L ||
    nrow(fsetdiff(expected, actual)) || nrow(fsetdiff(actual, expected))) {
  stop("DIAGNOSTIC_TEN_COHORT_WINDOWS_REQUIRED", call. = FALSE)
}
outputs <- c(pdf_path, rows_path, receipt_path)
if (any(file.exists(outputs))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(outputs)))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
input_name <- if (kind == "proportional") "agreement_proportional_bias.tsv" else "agreement_heteroscedasticity.tsv"
records <- lapply(seq_len(nrow(manifest)), function(i) {
  path <- file.path(manifest$diagnostic_dir[[i]], input_name)
  if (!file.exists(path)) stop("DIAGNOSTIC_SOURCE_MISSING", call. = FALSE)
  row <- fread(path)
  needed <- if (kind == "proportional") {
    c("cohort", "window_minutes", "slope", "ci_lower", "ci_upper", "model_estimated_bias_at_90")
  } else {
    c("cohort", "window_minutes", "variance_function_parameter", "ci_lower", "ci_upper",
      "residual_sd_at_90", "likelihood_ratio", "p_value")
  }
  if (nrow(row) != 1L || !all(needed %in% names(row))) stop("DIAGNOSTIC_SOURCE_CONTRACT_INVALID", call. = FALSE)
  if (row$cohort[[1L]] != manifest$cohort[[i]] ||
      row$window_minutes[[1L]] != manifest$window_minutes[[i]]) {
    stop("DIAGNOSTIC_SOURCE_IDENTITY_MISMATCH", call. = FALSE)
  }
  values <- unlist(row[, ..needed][, -c("cohort", "window_minutes"), with = FALSE], use.names = FALSE)
  if (any(!is.finite(as.numeric(values)))) stop("DIAGNOSTIC_NONFINITE_VALUE", call. = FALSE)
  if (kind == "proportional") {
    if (row$ci_lower > row$slope || row$slope > row$ci_upper) stop("PROPORTIONAL_CI_INVALID", call. = FALSE)
  } else {
    if (row$ci_lower > row$variance_function_parameter ||
        row$variance_function_parameter > row$ci_upper || row$residual_sd_at_90 <= 0 ||
        row$likelihood_ratio < 0 || row$p_value < 0 || row$p_value > 1) {
      stop("HETEROSCEDASTIC_VALUES_INVALID", call. = FALSE)
    }
  }
  row
})
data <- rbindlist(records, use.names = TRUE, fill = TRUE)
setkey(data, cohort, window_minutes)
cohort_names <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
                  SICDB = "SICdb", Lianyungang = "Lianyungang")
fmt_ci <- function(est, lo, hi, digits) {
  pattern <- paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)")
  sprintf(pattern, est, lo, hi)
}
rows <- rbindlist(lapply(c(60L, 5L), function(window_value) {
  rbindlist(lapply(cohorts, function(cohort_value) {
    z <- data[cohort == cohort_value & window_minutes == window_value]
    if (nrow(z) != 1L) stop("DIAGNOSTIC_LOOKUP_NOT_UNIQUE", call. = FALSE)
    if (kind == "proportional") {
      data.table(cohort = unname(cohort_names[[cohort_value]]), window = paste0(window_value, " min"),
                 slope_with_95ci = fmt_ci(z$slope, z$ci_lower, z$ci_upper, 2L),
                 model_estimated_bias_at_90 = sprintf("%.2f", z$model_estimated_bias_at_90))
    } else {
      data.table(cohort = unname(cohort_names[[cohort_value]]), window = paste0(window_value, " min"),
                 variance_parameter_with_95ci = fmt_ci(z$variance_function_parameter,
                                                        z$ci_lower, z$ci_upper, 3L),
                 residual_sd_at_90 = sprintf("%.2f", z$residual_sd_at_90),
                 likelihood_ratio = sprintf("%.1f", z$likelihood_ratio),
                 p_value = if (z$p_value < 0.001) "<.001" else sprintf("%.3f", z$p_value))
    }
  }))
}))
if (nrow(rows) != 10L) stop(paste("DIAGNOSTIC_RENDER_ROW_COUNT_INVALID", nrow(rows)), call. = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
grid.newpage()
title <- if (kind == "proportional") {
  "Supplementary Table 14. Continuous Proportional Bias Diagnostics"
} else {
  "Supplementary Table 15. Heteroscedasticity Diagnostics"
}
grid.text(title, x = unit(0.4, "in"), y = unit(7.72, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 11, fontface = "bold"))
if (kind == "proportional") {
  xs <- c(0.5, 2.9, 4.4, 7.8)
  headers <- c("Cohort", "Window", "Slope per percentage point (95% CI)",
               "Model-estimated Bias at mean saturation 90%")
} else {
  xs <- c(0.42, 2.25, 3.55, 6.3, 8.8, 10.15)
  headers <- c("Cohort", "Window", "Variance parameter (95% CI)",
               "Residual SD at 90%", "Likelihood ratio", "P value")
}
grid.lines(x = unit(c(0.4, 11.3), "in"), y = unit(c(6.99, 6.99), "in"), gp = gpar(lwd = 1.2))
for (j in seq_along(headers)) {
  grid.text(headers[[j]], x = unit(xs[[j]], "in"), y = unit(6.64, "in"),
            just = c("left", "centre"), gp = gpar(fontsize = 7.8, fontface = "bold"))
}
grid.lines(x = unit(c(0.4, 11.3), "in"), y = unit(c(6.34, 6.34), "in"), gp = gpar(lwd = 0.7))
ys <- seq(6.03, by = -0.4, length.out = 10L)
for (i in seq_len(nrow(rows))) {
  if (i == 6L) grid.lines(x = unit(c(0.4, 11.3), "in"),
                          y = unit(c(ys[[5L]] - 0.2, ys[[5L]] - 0.2), "in"),
                          gp = gpar(lwd = 0.6, col = "grey55"))
  for (j in seq_len(ncol(rows))) {
    grid.text(as.character(rows[[j]][[i]]), x = unit(xs[[j]], "in"), y = unit(ys[[i]], "in"),
              just = c("left", "centre"), gp = gpar(fontsize = 8.0))
  }
}
grid.lines(x = unit(c(0.4, 11.3), "in"), y = unit(c(2.18, 2.18), "in"), gp = gpar(lwd = 1.1))
note <- if (kind == "proportional") {
  "Slope is change in SpO2-SaO2 Bias per percentage point of paired mean saturation."
} else {
  "The likelihood ratio compares constant and saturation-varying residual-variance models."
}
grid.text(note, x = unit(0.42, "in"), y = unit(1.87, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 7.5))
dev.off()
fwrite(rows, rows_path, sep = "\t")
fwrite(data.table(table_id = if (kind == "proportional") "Supplementary Table 14" else "Supplementary Table 15",
                  input_n = nrow(data), row_n = nrow(rows), cohort_n = 5L, analysis_n = 10L),
       receipt_path, sep = "\t")
cat("AGREEMENT_DIAGNOSTIC_TABLE_RENDER_PASS kind=", kind, "\n", sep = "")
