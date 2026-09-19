#!/usr/bin/env Rscript
# Candidate five-cohort ICU-day figure from validated threshold aggregates.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: MANIFEST.tsv WINDOW(60|5) NEW_OUTPUT.pdf NEW_RECEIPT.tsv", call. = FALSE)
manifest <- fread(args[[1L]])
window <- suppressWarnings(as.integer(args[[2L]]))
pdf_path <- args[[3L]]
receipt_path <- args[[4L]]
if (is.na(window) || !(window %in% c(60L, 5L))) stop("WINDOW_INVALID", call. = FALSE)
if (file.exists(pdf_path) || file.exists(receipt_path)) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (!dir.exists(dirname(pdf_path)) || !dir.exists(dirname(receipt_path))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
if (!identical(names(manifest), c("cohort", "day_tsv"))) stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (!setequal(manifest$cohort, cohorts) || nrow(manifest) != 5L || anyDuplicated(manifest$cohort)) {
  stop("FIVE_COHORT_SET_INVALID", call. = FALSE)
}
if (any(!file.exists(manifest$day_tsv))) stop("MANIFEST_INPUT_MISSING", call. = FALSE)
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
data <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  x <- fread(manifest$day_tsv[[i]])
  required <- c("cohort", "window_minutes", "relative_icu_day", "spo2_stratum",
                "numerator_n", "denominator_n", "patient_n", "conditional_proportion",
                "ci_lower", "ci_upper", "cell_status")
  if (!all(required %in% names(x))) stop("DAY_INPUT_COLUMNS_MISSING", call. = FALSE)
  if (nrow(x) != 28L || any(x$cohort != manifest$cohort[[i]]) || any(x$window_minutes != window)) {
    stop(paste("COHORT_WINDOW_OR_CELL_COUNT_MISMATCH", manifest$cohort[[i]]), call. = FALSE)
  }
  x
}), use.names = TRUE, fill = TRUE)
expected_grid <- CJ(cohort = cohorts, relative_icu_day = 1:7, spo2_stratum = strata)
actual_grid <- unique(data[, .(cohort, relative_icu_day, spo2_stratum)])
if (nrow(data) != 140L || nrow(actual_grid) != 140L ||
    nrow(fsetdiff(expected_grid, actual_grid)) != 0L ||
    nrow(fsetdiff(actual_grid, expected_grid)) != 0L) {
  stop("DAY_CELL_GRID_INVALID", call. = FALSE)
}
if (anyNA(data[, .(numerator_n, denominator_n, patient_n, cell_status)]) ||
    any(!data$cell_status %in% c("DISPLAYABLE", "SPARSE", "NO_DENOMINATOR")) ||
    any(data$numerator_n < 0 | data$denominator_n < 0 | data$patient_n < 0 |
        data$numerator_n > data$denominator_n | data$patient_n > data$denominator_n)) {
  stop("DAY_COUNT_OR_STATUS_INVALID", call. = FALSE)
}
if (data[cell_status == "DISPLAYABLE", any(denominator_n < 20 | patient_n < 10)]) {
  stop("DISPLAYABLE_CELL_IS_SPARSE", call. = FALSE)
}
if (data[cell_status == "DISPLAYABLE", any(is.na(conditional_proportion) | is.na(ci_lower) | is.na(ci_upper) |
                                           conditional_proportion < 0 | conditional_proportion > 1 |
                                           ci_lower < 0 | ci_upper > 1 |
                                           ci_lower > conditional_proportion | ci_upper < conditional_proportion)]) {
  stop("DAY_ESTIMATE_OR_CI_INVALID", call. = FALSE)
}
if (data[cell_status != "DISPLAYABLE", any(denominator_n >= 20 & patient_n >= 10)]) {
  stop("SPARSE_CELL_STATUS_INVALID", call. = FALSE)
}
if (data[cell_status == "NO_DENOMINATOR",
         any(denominator_n != 0L | patient_n != 0L | numerator_n != 0L |
             !is.na(conditional_proportion) | !is.na(ci_lower) | !is.na(ci_upper))]) {
  stop("NO_DENOMINATOR_CELL_INVALID", call. = FALSE)
}

offsets <- setNames(c(-0.15, -0.05, 0.05, 0.15), strata)
line_types <- setNames(c("solid", "22", "42", "13"), strata)
shapes <- setNames(c(16, 17, 15, 18), strata)
colors <- setNames(c("#0072B2", "#D55E00", "#009E73", "#7A5195"), strata)
labels <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
data[, `:=`(stratum_factor = factor(spo2_stratum, levels = strata),
             plot_x = relative_icu_day + offsets[spo2_stratum],
             estimate_pct = fifelse(cell_status == "DISPLAYABLE", 100 * conditional_proportion, NA_real_),
             lower_pct = fifelse(cell_status == "DISPLAYABLE", 100 * ci_lower, NA_real_),
             upper_pct = fifelse(cell_status == "DISPLAYABLE", 100 * ci_upper, NA_real_))]
display_data <- data[cell_status == "DISPLAYABLE"]
line_data <- copy(display_data)
setorder(line_data, cohort, spo2_stratum, relative_icu_day)
line_data[, segment := cumsum(c(TRUE, diff(relative_icu_day) != 1L)), by = .(cohort, spo2_stratum)]
line_data[, segment_n := .N, by = .(cohort, spo2_stratum, segment)]
line_data <- line_data[segment_n >= 2L]
lower_observed <- data[cell_status == "DISPLAYABLE" & spo2_stratum %chin% strata[3:4], max(upper_pct, na.rm = TRUE)]
lower_limit <- if (is.finite(lower_observed) && lower_observed <= 10) 10 else
  min(100, max(15, ceiling(lower_observed / 5) * 5))
lower_breaks <- seq(0, lower_limit, by = if (lower_limit <= 10) 2 else 5)

base_theme <- theme_classic(base_size = 10) + theme(
  axis.line = element_line(linewidth = 0.35, color = "black"),
  axis.text = element_text(size = 8.5, color = "black"),
  panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
  panel.grid.minor = element_blank(), plot.title = element_text(size = 10.5, face = "bold"),
  legend.position = "none")
subpanel <- function(cohort_name, chosen, y_limit, y_breaks, title, show_x) {
  p <- ggplot(display_data[cohort == cohort_name & spo2_stratum %chin% chosen],
              aes(x = plot_x, y = estimate_pct, color = stratum_factor,
                  linetype = stratum_factor, shape = stratum_factor)) +
    geom_errorbar(aes(ymin = lower_pct, ymax = upper_pct), width = 0.045,
                  linewidth = 0.28, alpha = 0.42) +
    geom_line(data = line_data[cohort == cohort_name & spo2_stratum %chin% chosen],
              aes(group = interaction(stratum_factor, segment)), linewidth = 0.58) +
    geom_point(size = 2) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_linetype_manual(values = line_types, drop = FALSE) +
    scale_shape_manual(values = shapes, drop = FALSE) +
    scale_x_continuous(breaks = 1:7, limits = c(0.72, 7.28), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = y_breaks, labels = function(z) paste0(z, "%"),
                       limits = c(0, y_limit), expand = expansion(mult = c(0, 0.015))) +
    labs(title = if (title) unname(labels[[cohort_name]]) else NULL, x = NULL, y = NULL) + base_theme
  if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  p
}
panels <- lapply(cohorts, function(cohort) {
  upper <- subpanel(cohort, strata[1:2], 100, seq(20, 100, 20), TRUE, FALSE)
  lower <- subpanel(cohort, strata[3:4], lower_limit, lower_breaks, FALSE, TRUE)
  upper / lower + plot_layout(heights = c(2, 1))
})
legend_labels <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
legend_data <- data.frame(stratum = factor(strata, levels = strata),
                          label = legend_labels, y = c(4.2, 3.5, 2.4, 1.7))
legend_panel <- ggplot(legend_data) +
  geom_segment(aes(x = 0.08, xend = 0.30, y = y, yend = y, linetype = stratum, color = stratum), linewidth = 0.65) +
  geom_point(aes(x = 0.19, y = y, shape = stratum, color = stratum), size = 2.1) +
  geom_text(aes(x = 0.36, y = y, label = label), hjust = 0, size = 3.25) +
  annotate("text", x = 0.08, y = 4.85, label = "SpO2 stratum", hjust = 0, fontface = "bold", size = 3.6) +
  annotate("text", x = 0.08, y = 4.55, label = "Upper panel: 0%-100%", hjust = 0, size = 2.8) +
  annotate("text", x = 0.08, y = 2.75, label = paste0("Lower panel: 0%-", lower_limit, "%"), hjust = 0, size = 2.8) +
  annotate("text", x = 0.08, y = 0.8, label = "Sparse or empty cells\nappear as gaps.",
           hjust = 0, vjust = 1, size = 2.8) +
  scale_color_manual(values = colors) + scale_linetype_manual(values = line_types) +
  scale_shape_manual(values = shapes) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5.1), clip = "off") +
  theme_void() + theme(legend.position = "none")
combined <- wrap_plots(c(panels, list(legend_panel)), ncol = 3, byrow = TRUE)
final <- ggdraw(combined) +
  draw_label("ICU day", x = 0.50, y = 0.012, hjust = 0.5, vjust = 0, size = 10.5) +
  draw_label(expression("Pairs with " * SaO[2] * " <88% (%)"), x = 0.017, y = 0.50,
             angle = 90, hjust = 0.5, vjust = 0.5, size = 10.5)
ggsave(pdf_path, final, width = 12, height = 9, units = "in", device = cairo_pdf, bg = "white")
receipt <- data[, .(cell_n = .N, displayable_n = sum(cell_status == "DISPLAYABLE"),
                    sparse_n = sum(cell_status != "DISPLAYABLE")), by = cohort]
receipt[, `:=`(window_minutes = window, lower_panel_max_pct = lower_limit)]
fwrite(receipt, receipt_path, sep = "\t")
cat("ICU_DAY_CANDIDATE_RENDER_PASS window=", window, " cells=", nrow(data), "\n", sep = "")
