#!/usr/bin/env Rscript
# Technical candidate for Supplementary Table 22, all model terms in both windows.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: LMM60_MANIFEST.tsv LMM5_MANIFEST.tsv NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
outputs <- args[3:5]
if (any(file.exists(outputs))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(outputs)))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
panels <- list()
for (window_index in seq_along(c(60L, 5L))) {
  window_value <- c(60L, 5L)[[window_index]]
  m <- fread(args[[window_index]])
  if (!identical(names(m), c("cohort", "lmm_tsv")) || nrow(m) != 5L ||
      !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort)) {
    stop("LMM_BALANCE_MANIFEST_INVALID", call. = FALSE)
  }
  for (cohort_value in cohorts) {
    primary_path <- m[cohort == cohort_value]$lmm_tsv
    balanced_path <- file.path(dirname(primary_path), "lmm_patient_balanced_pooled.tsv")
    if (length(primary_path) != 1L || !file.exists(primary_path) || !file.exists(balanced_path)) {
      stop("LMM_BALANCE_SOURCE_MISSING", call. = FALSE)
    }
    p <- fread(primary_path)
    b <- fread(balanced_path)
    needed <- c("term", "estimate", "ci_lower", "ci_upper", "reporting_label", "reporting_multiplier")
    if (!all(needed %in% names(p)) || !all(needed %in% names(b)) ||
        nrow(p) < 1L || nrow(p) != nrow(b) || anyDuplicated(p$term) ||
        anyDuplicated(b$term) || !setequal(p$term, b$term)) {
      stop("LMM_BALANCE_TERM_CONTRACT_INVALID", call. = FALSE)
    }
    b <- b[match(p$term, b$term)]
    if (anyNA(p[, ..needed]) || anyNA(b[, ..needed]) ||
        any(!is.finite(as.matrix(p[, .(estimate, ci_lower, ci_upper, reporting_multiplier)]))) ||
        any(!is.finite(as.matrix(b[, .(estimate, ci_lower, ci_upper, reporting_multiplier)]))) ||
        any(p$reporting_multiplier <= 0 | b$reporting_multiplier <= 0) ||
        any(p$reporting_multiplier != b$reporting_multiplier | p$reporting_label != b$reporting_label) ||
        any(p$ci_lower > p$estimate | p$estimate > p$ci_upper |
            b$ci_lower > b$estimate | b$estimate > b$ci_upper)) {
      stop("LMM_BALANCE_VALUE_CONTRACT_INVALID", call. = FALSE)
    }
    fmt <- function(est, lo, hi, multiplier) {
      sprintf("%+.2f (%+.2f to %+.2f)", est * multiplier, lo * multiplier, hi * multiplier)
    }
    panel <- data.table(cohort = unname(pretty[[cohort_value]]), window_minutes = window_value,
      term = p$term, label = p$reporting_label,
      primary_estimate = p$estimate * p$reporting_multiplier,
      balanced_estimate = b$estimate * b$reporting_multiplier,
      primary_with_95ci = fmt(p$estimate, p$ci_lower, p$ci_upper, p$reporting_multiplier),
      balanced_with_95ci = fmt(b$estimate, b$ci_lower, b$ci_upper, b$reporting_multiplier),
      same_direction = ifelse(sign(p$estimate) == sign(b$estimate), "Yes", "No"))
    panels[[paste(cohort_value, window_value)]] <- panel
  }
}
rows <- rbindlist(panels)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(outputs[[1L]], width = 11.69, height = 8.27, family = "sans", bg = "white")
for (window_value in c(60L, 5L)) for (cohort_value in cohorts) {
  panel <- panels[[paste(cohort_value, window_value)]]
  if (nrow(panel) > 25L) stop("LMM_BALANCE_PANEL_TOO_LONG_FOR_TECHNICAL_PAGE", call. = FALSE)
  grid.newpage()
  grid.text("Supplementary Table 22. Primary and Patient-Balanced LMM Estimates",
    x = unit(0.42, "in"), y = unit(7.72, "in"), just = c("left", "centre"),
    gp = gpar(fontsize = 10.5, fontface = "bold"))
  grid.text(paste(unname(pretty[[cohort_value]]), "-", window_value, "min"),
    x = unit(0.42, "in"), y = unit(7.31, "in"), just = c("left", "centre"),
    gp = gpar(fontsize = 9, fontface = "bold"))
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(7.02, 7.02), "in"), gp = gpar(lwd = 1.0))
  xs <- c(0.46, 5.35, 7.95, 10.7)
  headers <- c("Variable, comparison, or unit", "Primary coefficient (95% CI)",
               "Patient-balanced coefficient (95% CI)", "Same direction")
  for (j in seq_along(headers)) grid.text(headers[[j]], x = unit(xs[[j]], "in"),
    y = unit(6.74, "in"), just = c("left", "centre"),
    gp = gpar(fontsize = 7.4, fontface = "bold"))
  grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.47, 6.47), "in"), gp = gpar(lwd = 0.55))
  ys <- seq(6.24, by = -0.22, length.out = nrow(panel))
  for (i in seq_len(nrow(panel))) {
    values <- unlist(panel[i, .(label, primary_with_95ci, balanced_with_95ci, same_direction)],
                     use.names = FALSE)
    for (j in seq_along(values)) grid.text(values[[j]], x = unit(xs[[j]], "in"),
      y = unit(ys[[i]], "in"), just = c("left", "centre"), gp = gpar(fontsize = 7.2))
  }
  grid.lines(x = unit(c(0.42, 11.27), "in"),
    y = unit(c(min(ys) - 0.18, min(ys) - 0.18), "in"), gp = gpar(lwd = 1.0))
}
dev.off()
fwrite(rows, outputs[[2L]], sep = "\t")
fwrite(data.table(table_id = "Supplementary Table 22", cohort_n = 5L, analysis_n = 10L,
                  panel_n = length(panels), row_n = nrow(rows)), outputs[[3L]], sep = "\t")
cat("LMM_BALANCE_TABLE_CANDIDATE_RENDER_PASS rows=", nrow(rows), "\n", sep = "")
