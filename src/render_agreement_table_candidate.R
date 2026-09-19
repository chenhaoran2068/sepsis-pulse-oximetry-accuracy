#!/usr/bin/env Rscript
# Five-cohort agreement table; numerical inputs come only from calculation outputs.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: MANIFEST.tsv WINDOW(60|5) NEW_OUTPUT.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
}
manifest <- fread(args[[1L]])
window <- suppressWarnings(as.integer(args[[2L]]))
pdf_path <- args[[3L]]
rows_path <- args[[4L]]
receipt_path <- args[[5L]]
if (is.na(window) || !window %in% c(60L, 5L)) stop("WINDOW_INVALID", call. = FALSE)
if (any(file.exists(c(pdf_path, rows_path, receipt_path)))) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(c(pdf_path, rows_path, receipt_path))))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
if (!all(c("cohort", "agreement_tsv") %in% names(manifest))) stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (nrow(manifest) != 5L || !setequal(manifest$cohort, cohorts) || anyDuplicated(manifest$cohort)) {
  stop("FIVE_COHORT_SET_INVALID", call. = FALSE)
}
if (any(!file.exists(manifest$agreement_tsv))) stop("MANIFEST_INPUT_MISSING", call. = FALSE)
required <- c("cohort", "window_minutes", "pair_n", "stay_n", "mean_bias",
              "mean_bias_ci_lower", "mean_bias_ci_upper", "lower_loa", "lower_loa_ci_lower",
              "lower_loa_ci_upper", "upper_loa", "upper_loa_ci_lower", "upper_loa_ci_upper",
              "arms", "arms_ci_lower", "arms_ci_upper")
data <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  x <- fread(manifest$agreement_tsv[[i]])
  if (!all(required %in% names(x))) stop("AGREEMENT_COLUMNS_MISSING", call. = FALSE)
  if (nrow(x) != 1L || x$cohort[[1L]] != manifest$cohort[[i]] ||
      x$window_minutes[[1L]] != window) {
    stop(paste("COHORT_WINDOW_OR_ROW_MISMATCH", manifest$cohort[[i]]), call. = FALSE)
  }
  x
}), use.names = TRUE, fill = TRUE)
data <- data[match(cohorts, cohort)]
if (anyNA(data[, ..required]) || any(data$pair_n < 1 | data$stay_n < 1 |
                                    data$stay_n > data$pair_n | data$arms < 0)) {
  stop("AGREEMENT_VALUE_INVALID", call. = FALSE)
}
for (prefix in c("mean_bias", "lower_loa", "upper_loa", "arms")) {
  center <- data[[prefix]]
  low <- data[[paste0(prefix, "_ci_lower")]]
  high <- data[[paste0(prefix, "_ci_upper")]]
  if (any(!is.finite(center) | !is.finite(low) | !is.finite(high) |
          low > center | center > high)) stop(paste("AGREEMENT_CI_INVALID", prefix), call. = FALSE)
}
if (any(data$lower_loa >= data$mean_bias | data$mean_bias >= data$upper_loa)) {
  stop("LOA_ORDER_INVALID", call. = FALSE)
}
pretty_name <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
                 SICDB = "SICdb", Lianyungang = "Lianyungang")
fmt_ci <- function(z, stem) sprintf("%.2f (%.2f to %.2f)", z[[stem]],
                                  z[[paste0(stem, "_ci_lower")]],
                                  z[[paste0(stem, "_ci_upper")]])
rows <- data.table(
  cohort = unname(pretty_name[data$cohort]),
  pairs_and_stays = sprintf("%s (%s)", format(data$pair_n, big.mark = ",", scientific = FALSE,
                                              trim = TRUE),
                            format(data$stay_n, big.mark = ",", scientific = FALSE, trim = TRUE)),
  mean_bias_ci = fmt_ci(data, "mean_bias"),
  lower_loa_ci = fmt_ci(data, "lower_loa"),
  upper_loa_ci = fmt_ci(data, "upper_loa"),
  arms_ci = fmt_ci(data, "arms")
)
title <- if (window == 60L) {
  "Table 1. Cohort-Specific Agreement Between SpO2 and SaO2 in the 60-Minute Primary Analysis"
} else {
  "Supplementary Table 11. Cohort-Specific Agreement Between SpO2 and SaO2 in the Independent 5-Minute Sensitivity Analysis"
}
headers <- c("Cohort", "Pairs, n (ICU stays, n)", "Mean Bias (95% CI)",
             "Lower LoA (95% CI)", "Upper LoA (95% CI)", "A_rms (95% CI)")
cell_x <- c(0.48, 1.95, 3.77, 5.69, 7.70, 9.68)
header_y <- 6.54
row_y <- seq(6.02, 4.02, by = -0.50)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
grid.newpage()
grid.text(title, x = unit(0.42, "in"), y = unit(7.70, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 11, fontface = "bold"))
grid.lines(x = unit(c(0.42, 11.28), "in"), y = unit(c(6.81, 6.81), "in"),
           gp = gpar(lwd = 1.3))
for (j in seq_along(headers)) {
  grid.text(headers[[j]], x = unit(cell_x[[j]], "in"), y = unit(header_y, "in"),
            just = c("left", "centre"), gp = gpar(fontsize = if (j == 2L) 8.0 else 8.8,
                                                   fontface = "bold"))
}
grid.lines(x = unit(c(0.42, 11.28), "in"), y = unit(c(6.29, 6.29), "in"), gp = gpar(lwd = 0.7))
for (i in seq_len(nrow(rows))) {
  for (j in seq_len(ncol(rows))) {
    grid.text(as.character(rows[[j]][[i]]), x = unit(cell_x[[j]], "in"),
              y = unit(row_y[[i]], "in"), just = c("left", "centre"),
              gp = gpar(fontsize = if (j == 1L) 8.7 else 8.1))
  }
}
grid.lines(x = unit(c(0.42, 11.28), "in"), y = unit(c(3.69, 3.69), "in"),
           gp = gpar(lwd = 1.1))
notes <- c(
  "Bias is SpO2 minus SaO2. Positive values indicate that SpO2 exceeds SaO2.",
  "The second column shows eligible pairs and contributing ICU stays. LoA uses the cohort-specific mixed-effects variance structure.",
  "Abbreviations: A_rms, root mean square accuracy. CI, confidence interval. ICU, intensive care unit. LoA, limits of agreement.",
  "AmsterdamUMCdb, Amsterdam University Medical Centers Database; eICU, eICU Collaborative Research Database;",
  "MIMIC, Medical Information Mart for Intensive Care; SICdb, Salzburg Intensive Care database;",
  "SaO2, arterial oxygen saturation. SpO2, peripheral oxygen saturation."
)
for (i in seq_along(notes)) {
  grid.text(notes[[i]], x = unit(0.43, "in"), y = unit(3.27 - (i - 1L) * 0.31, "in"),
            just = c("left", "centre"), gp = gpar(fontsize = 7.5))
}
dev.off()
fwrite(rows, rows_path, sep = "\t")
receipt <- data.table(window_minutes = window, table_id = if (window == 60L) "Main Table 1" else "Supplementary Table 11",
                      cohort_n = 5L, row_n = nrow(rows), pair_n_total = sum(data$pair_n))
fwrite(receipt, receipt_path, sep = "\t")
cat("AGREEMENT_TABLE_CANDIDATE_RENDER_PASS window=", window, " cohorts=5\n", sep = "")
