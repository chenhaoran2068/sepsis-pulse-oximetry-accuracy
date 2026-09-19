#!/usr/bin/env Rscript
# Main Table 2 from the ten independent cohort/window threshold aggregates.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: MANIFEST.tsv NEW_OUTPUT.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
manifest <- fread(args[[1L]])
pdf_path <- args[[2L]]
rows_path <- args[[3L]]
receipt_path <- args[[4L]]
if (!identical(names(manifest), c("cohort", "window_minutes", "low_sao2_tsv"))) {
  stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
}
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
expected <- CJ(cohort = cohorts, window_minutes = c(5L, 60L))
actual <- unique(manifest[, .(cohort, window_minutes)])
if (nrow(manifest) != 10L || nrow(actual) != 10L ||
    nrow(fsetdiff(expected, actual)) || nrow(fsetdiff(actual, expected))) {
  stop("TEN_COHORT_WINDOW_SET_INVALID", call. = FALSE)
}
if (any(!file.exists(manifest$low_sao2_tsv))) stop("MANIFEST_INPUT_MISSING", call. = FALSE)
if (any(file.exists(c(pdf_path, rows_path, receipt_path)))) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(c(pdf_path, rows_path, receipt_path))))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
needed <- c("threshold_spo2_ge", "accepted_pair_n", "pair_denominator_sao2_lt88_n",
            "pair_numerator_n", "pair_conditional_proportion", "pair_ci_lower", "pair_ci_upper",
            "cell_status")
data <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  x <- fread(manifest$low_sao2_tsv[[i]])
  if (!all(needed %in% names(x))) stop("THRESHOLD_COLUMNS_MISSING", call. = FALSE)
  if (nrow(x) != 2L || !setequal(x$threshold_spo2_ge, c(88, 92))) {
    stop("THRESHOLD_TWO_ROWS_INVALID", call. = FALSE)
  }
  x[, `:=`(cohort = manifest$cohort[[i]], window_minutes = manifest$window_minutes[[i]])]
  x
}), use.names = TRUE, fill = TRUE)
count_fields <- c("threshold_spo2_ge", "accepted_pair_n", "pair_denominator_sao2_lt88_n",
                  "pair_numerator_n", "cell_status")
if (anyNA(data[, ..count_fields]) ||
    any(data$accepted_pair_n < 1 | data$pair_denominator_sao2_lt88_n < 1 |
        data$pair_denominator_sao2_lt88_n > data$accepted_pair_n |
        data$pair_numerator_n < 0 | data$pair_numerator_n > data$pair_denominator_sao2_lt88_n) ||
    any(!data$cell_status %in% c("ESTIMATED", "SPARSE"))) {
  stop("THRESHOLD_COUNTS_OR_STATUS_INVALID", call. = FALSE)
}
if (data[cell_status == "ESTIMATED", any(!is.finite(pair_conditional_proportion) |
        !is.finite(pair_ci_lower) | !is.finite(pair_ci_upper) |
        abs(pair_conditional_proportion - pair_numerator_n / pair_denominator_sao2_lt88_n) > 1e-10 |
        pair_ci_lower < 0 | pair_ci_upper > 1 |
        pair_ci_lower > pair_conditional_proportion | pair_ci_upper < pair_conditional_proportion)]) {
  stop("THRESHOLD_ESTIMATE_OR_CI_INVALID", call. = FALSE)
}
by_group <- data[, .(denom_unique = uniqueN(pair_denominator_sao2_lt88_n),
                     n88 = pair_numerator_n[threshold_spo2_ge == 88],
                     n92 = pair_numerator_n[threshold_spo2_ge == 92]),
                 by = .(cohort, window_minutes)]
if (any(by_group$denom_unique != 1L | by_group$n92 > by_group$n88)) {
  stop("THRESHOLD_NESTING_INVALID", call. = FALSE)
}
data <- data[order(match(cohort, cohorts), -window_minutes, threshold_spo2_ge)]
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_episode <- function(x) {
  nd <- paste0(fmt_n(x$pair_numerator_n), "/", fmt_n(x$pair_denominator_sao2_lt88_n))
  if (x$cell_status == "SPARSE") nd else sprintf("%s (%.1f)", nd, 100 * x$pair_conditional_proportion)
}
fmt_ci <- function(x) {
  if (x$cell_status == "SPARSE") "-" else sprintf("%.1f to %.1f", 100 * x$pair_ci_lower,
                                                   100 * x$pair_ci_upper)
}
labels <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
out <- rbindlist(lapply(seq_along(cohorts), function(i) {
  rows <- lapply(c(60L, 5L), function(w) {
    x <- data[cohort == cohorts[[i]] & window_minutes == w]
    a <- x[threshold_spo2_ge == 88]
    b <- x[threshold_spo2_ge == 92]
    data.table(cohort = if (w == 60L) unname(labels[[cohorts[[i]]]]) else "",
               window = paste0(w, " min"),
               sao2_lt88_n = fmt_n(a$pair_denominator_sao2_lt88_n),
               occult_n_N_pct = fmt_episode(a), occult_ci_pct = fmt_ci(a),
               severe_n_N_pct = fmt_episode(b), severe_ci_pct = fmt_ci(b))
  })
  rbindlist(rows)
}))
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
grid.newpage()
grid.text("Table 2. Observed Occult and Severe Occult Hypoxemia by Cohort and Pairing Window",
          x = unit(0.42, "in"), y = unit(7.73, "in"), just = c("left", "centre"),
          gp = gpar(fontsize = 11, fontface = "bold"))
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.88, 6.88), "in"), gp = gpar(lwd = 1.2))
col_x <- c(0.48, 1.82, 2.75, 4.14, 6.08, 7.47, 9.56)
heads <- c("Cohort", "Window", "SaO2 <88%, n", "Occult, n/N (%)", "95% CI, %",
           "Severe occult, n/N (%)", "95% CI, %")
for (j in seq_along(heads)) {
  grid.text(heads[[j]], x = unit(col_x[[j]], "in"), y = unit(6.58, "in"),
            just = c("left", "centre"), gp = gpar(fontsize = 8.1, fontface = "bold"))
}
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.30, 6.30), "in"), gp = gpar(lwd = 0.7))
row_y <- seq(6.02, 2.42, by = -0.40)
for (i in seq_len(nrow(out))) {
  if (i > 1L && i %% 2L == 1L) {
    yy <- (row_y[[i - 1L]] + row_y[[i]]) / 2
    grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(yy, yy), "in"),
               gp = gpar(lwd = 0.35, col = "grey60"))
  }
  for (j in seq_len(ncol(out))) {
    grid.text(as.character(out[[j]][[i]]), x = unit(col_x[[j]], "in"),
              y = unit(row_y[[i]], "in"), just = c("left", "centre"),
              gp = gpar(fontsize = 8.0))
  }
}
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(2.18, 2.18), "in"), gp = gpar(lwd = 1.1))
notes <- c(
  "The denominator is eligible pairs with SaO2 <88%. Occult and severe occult hypoxemia require SpO2 at least 88% and 92%, respectively.",
  "Confidence intervals use 2,000 patient-clustered bootstrap resamples. The 60-minute and 5-minute pairs were constructed independently.",
  "Abbreviations: AmsterdamUMCdb, Amsterdam University Medical Centers Database. eICU, eICU Collaborative Research Database.",
  "MIMIC, Medical Information Mart for Intensive Care. SICdb, Salzburg Intensive Care database. CI, confidence interval.",
  "SaO2, arterial oxygen saturation. SpO2, peripheral oxygen saturation."
)
for (i in seq_along(notes)) {
  grid.text(notes[[i]], x = unit(0.43, "in"), y = unit(1.88 - (i - 1L) * 0.26, "in"),
            just = c("left", "centre"), gp = gpar(fontsize = 7.1))
}
dev.off()
fwrite(out, rows_path, sep = "\t")
fwrite(data.table(table_id = "Main Table 2", cohort_n = 5L, analysis_n = 10L,
                  input_row_n = nrow(data), display_row_n = nrow(out)), receipt_path, sep = "\t")
cat("THRESHOLD_TABLE_CANDIDATE_RENDER_PASS cohorts=5 analyses=10\n")
