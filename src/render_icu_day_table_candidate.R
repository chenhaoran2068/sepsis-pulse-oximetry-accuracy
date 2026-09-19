#!/usr/bin/env Rscript
# Supplementary Tables 18 and 19 from complete ICU-day-by-SpO2 aggregates.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST.tsv WINDOW(60|5) NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
window <- suppressWarnings(as.integer(args[[2L]]))
pdf_path <- args[[3L]]
rows_path <- args[[4L]]
receipt_path <- args[[5L]]
if (is.na(window) || !window %in% c(60L, 5L)) stop("WINDOW_INVALID", call. = FALSE)
if (!identical(names(m), c("cohort", "day_tsv"))) stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
labels <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
if (nrow(m) != 5L || !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort)) {
  stop("FIVE_COHORT_SET_INVALID", call. = FALSE)
}
if (any(!file.exists(m$day_tsv))) stop("MANIFEST_INPUT_MISSING", call. = FALSE)
if (any(file.exists(c(pdf_path, rows_path, receipt_path)))) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(c(pdf_path, rows_path, receipt_path))))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
required <- c("cohort", "window_minutes", "relative_icu_day", "spo2_stratum", "denominator_n",
              "patient_n", "numerator_n", "conditional_proportion", "ci_lower", "ci_upper", "cell_status")
data <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
  x <- fread(m$day_tsv[[i]])
  if (!all(required %in% names(x))) stop("DAY_INPUT_COLUMNS_MISSING", call. = FALSE)
  if (nrow(x) != 28L || any(x$cohort != m$cohort[[i]]) || any(x$window_minutes != window)) {
    stop("COHORT_WINDOW_OR_CELL_COUNT_MISMATCH", call. = FALSE)
  }
  x
}), use.names = TRUE, fill = TRUE)
expected <- CJ(cohort = cohorts, relative_icu_day = 1:7, spo2_stratum = strata)
actual <- unique(data[, .(cohort, relative_icu_day, spo2_stratum)])
if (nrow(data) != 140L || nrow(actual) != 140L ||
    nrow(fsetdiff(expected, actual)) || nrow(fsetdiff(actual, expected))) {
  stop("DAY_CELL_GRID_INVALID", call. = FALSE)
}
if (anyNA(data[, .(denominator_n, patient_n, numerator_n, cell_status)]) ||
    any(!data$cell_status %in% c("DISPLAYABLE", "SPARSE", "NO_DENOMINATOR")) ||
    any(data$denominator_n < 0 | data$patient_n < 0 | data$numerator_n < 0 |
        data$patient_n > data$denominator_n | data$numerator_n > data$denominator_n)) {
  stop("DAY_COUNT_OR_STATUS_INVALID", call. = FALSE)
}
if (data[cell_status == "DISPLAYABLE", any(denominator_n < 20 | patient_n < 10 |
      is.na(conditional_proportion) | is.na(ci_lower) | is.na(ci_upper) |
      abs(conditional_proportion - numerator_n / denominator_n) > 1e-10 |
      ci_lower < 0 | ci_upper > 1 | ci_lower > conditional_proportion |
      ci_upper < conditional_proportion)]) {
  stop("DAY_ESTIMATE_OR_CI_INVALID", call. = FALSE)
}
if (data[cell_status == "SPARSE", any(denominator_n >= 20 & patient_n >= 10)]) {
  stop("SPARSE_CELL_STATUS_INVALID", call. = FALSE)
}
if (data[cell_status == "NO_DENOMINATOR",
         any(denominator_n != 0L | patient_n != 0L | numerator_n != 0L |
             !is.na(conditional_proportion) | !is.na(ci_lower) | !is.na(ci_upper))]) {
  stop("NO_DENOMINATOR_CELL_INVALID", call. = FALSE)
}
fmt_n <- function(z) format(z, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_cell <- function(z) {
  count <- paste0(fmt_n(z$numerator_n), "/", fmt_n(z$denominator_n))
  if (z$cell_status %in% c("SPARSE", "NO_DENOMINATOR")) return(count)
  sprintf("%s (%.1f) [%.1f to %.1f]", count, 100 * z$conditional_proportion,
          100 * z$ci_lower, 100 * z$ci_upper)
}
out <- rbindlist(lapply(seq_along(cohorts), function(i) {
  rbindlist(lapply(1:7, function(d) {
    x <- data[cohort == cohorts[[i]] & relative_icu_day == d]
    vals <- vapply(strata, function(s) fmt_cell(x[spo2_stratum == s]), character(1L))
    data.table(cohort = if (d == 1L) unname(labels[[cohorts[[i]]]]) else "",
               icu_day = as.character(d), spo2_70_to_lt88 = vals[[1L]],
               spo2_88_to_lt92 = vals[[2L]], spo2_92_to_le96 = vals[[3L]],
               spo2_gt96_to_100 = vals[[4L]])
  }))
}))
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
title <- paste0("Supplementary Table ", if (window == 60L) "18" else "19",
                ". ICU Day 1-7 Arterial Hypoxemia by Displayed SpO2 Stratum (",
                window, "-Minute Analysis)")
headers <- c("Cohort", "Day", "SpO2 70%-<88%", "SpO2 88%-<92%", "SpO2 92%-<=96%", "SpO2 >96%-100%")
xs <- c(0.42, 1.72, 2.20, 4.50, 6.72, 8.95)
draw_page <- function(first, last, page) {
  grid.newpage()
  grid.text(title, x = unit(0.42, "in"), y = unit(7.77, "in"), just = c("left", "centre"),
            gp = gpar(fontsize = 10.5, fontface = "bold"))
  grid.text(paste0("Page ", page), x = unit(11.22, "in"), y = unit(7.77, "in"),
            just = c("right", "centre"), gp = gpar(fontsize = 8))
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(7.35, 7.35), "in"), gp = gpar(lwd = 1.1))
  for (j in seq_along(headers)) {
    grid.text(headers[[j]], x = unit(xs[[j]], "in"), y = unit(7.12, "in"),
              just = c("left", "centre"), gp = gpar(fontsize = 7.4, fontface = "bold"))
  }
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.91, 6.91), "in"), gp = gpar(lwd = 0.6))
  ids <- first:last
  yy <- seq(6.69, by = -0.205, length.out = length(ids))
  for (k in seq_along(ids)) {
    i <- ids[[k]]
    if (k > 1L && (i - 1L) %% 7L == 0L) {
      line_y <- (yy[[k - 1L]] + yy[[k]]) / 2
      grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(line_y, line_y), "in"),
                 gp = gpar(lwd = 0.3, col = "grey60"))
    }
    for (j in seq_len(ncol(out))) {
      grid.text(as.character(out[[j]][[i]]), x = unit(xs[[j]], "in"), y = unit(yy[[k]], "in"),
                just = c("left", "centre"), gp = gpar(fontsize = 6.7))
    }
  }
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(2.13, 2.13), "in"), gp = gpar(lwd = 1.0))
  notes <- c("Each cell shows the count with SaO2 <88% / eligible pairs in that ICU-day-by-SpO2 stratum, percentage, and 95% CI.",
             "Cells with fewer than 20 eligible pairs or fewer than 10 patients show counts only.",
             "Abbreviations: CI, confidence interval. ICU, intensive care unit. SaO2, arterial oxygen saturation. SpO2, peripheral oxygen saturation.")
  for (i in seq_along(notes)) {
    grid.text(notes[[i]], x = unit(0.43, "in"), y = unit(1.79 - (i - 1L) * 0.30, "in"),
              just = c("left", "centre"), gp = gpar(fontsize = 7.0))
  }
}
draw_page(1L, 21L, 1L)
draw_page(22L, 35L, 2L)
dev.off()
fwrite(out, rows_path, sep = "\t")
fwrite(data.table(table_id = if (window == 60L) "Supplementary Table 18" else "Supplementary Table 19",
                  window_minutes = window, cohort_n = 5L, input_cell_n = nrow(data),
                  display_row_n = nrow(out), sparse_n = sum(data$cell_status != "DISPLAYABLE"), page_n = 2L),
       receipt_path, sep = "\t")
cat("ICU_DAY_TABLE_CANDIDATE_RENDER_PASS window=", window, " cells=140 rows=35\n", sep = "")
