#!/usr/bin/env Rscript
# Technical candidates for Supplementary Figures 7-11 from cohort-specific forest outputs.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: RF_DIR COHORT FIGURE_NUMBER NEW.pdf NEW_RECEIPT.tsv", call. = FALSE)
rf_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cohort <- args[[2L]]
figure_number <- suppressWarnings(as.integer(args[[3L]]))
pdf_path <- args[[4L]]
receipt_path <- args[[5L]]
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (!cohort %in% cohorts || is.na(figure_number) || figure_number != match(cohort, cohorts) + 6L) {
  stop("RF_SHAP_FIGURE_IDENTITY_INVALID", call. = FALSE)
}
if (any(file.exists(c(pdf_path, receipt_path)))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(c(pdf_path, receipt_path))))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
paths <- file.path(rf_dir, c("rf_shap_grouped.rds", "rf_shap_summary.tsv", "rf_shap_sample.tsv",
  "rf_shap_feature_values.tsv", "rf_diagnostics.tsv", "status.txt"))
if (any(!file.exists(paths))) stop("RF_SHAP_DISPLAY_SOURCE_MISSING", call. = FALSE)
names(paths) <- c("grouped", "summary", "sample", "values", "diagnostics", "status")
grouped <- readRDS(paths[["grouped"]])
summary <- fread(paths[["summary"]])
sample <- fread(paths[["sample"]])
values <- fread(paths[["values"]])
diagnostics <- fread(paths[["diagnostics"]])
status <- trimws(readLines(paths[["status"]], warn = FALSE)[[1L]])
if (!status %in% c("RF_PAPER_SETTINGS_COMPLETE_CANDIDATE_ONLY",
                   "RF_SMOKE_COMPLETE_NOT_PAPER_SETTINGS") ||
    nrow(diagnostics) != 1L || diagnostics$cohort[[1L]] != cohort ||
    !isTRUE(diagnostics$paired_mean_integerized[[1L]]) ||
    isTRUE(diagnostics$standalone_spo2_or_sao2[[1L]]) ||
    !is.matrix(grouped) || nrow(grouped) != nrow(sample) ||
    nrow(grouped) != nrow(values) || nrow(grouped) < 10L || nrow(grouped) > 1000L ||
    ncol(grouped) != nrow(summary) ||
    !identical(colnames(grouped), summary$feature_group) ||
    anyDuplicated(sample$analysis_patient_key) ||
    anyDuplicated(sample$pair_uid) ||
    !identical(as.character(sample$pair_uid), as.character(values$pair_uid)) ||
    any(!is.finite(grouped)) || any(!is.finite(summary$mean_abs_shap)) ||
    any(abs(colMeans(abs(grouped)) - summary$mean_abs_shap) > 1e-7) ||
    any(summary$feature_group %in% c("SpO2", "SaO2"))) {
  stop("RF_SHAP_DISPLAY_SOURCE_CONTRACT_INVALID", call. = FALSE)
}
feature_sources <- c("Paired mean saturation" = "paired_mean_saturation_c90_rf",
  "Age" = if ("age_value_per_10y" %in% names(values)) "age_value_per_10y" else "age_interval",
  "Sex" = "sex_class", "Charlson Comorbidity Index" = "cci_value",
  "SOFA" = "sofa_total", "PaO2" = "pao2_candidate_value",
  "PaCO2" = "paco2_candidate_value", "pH" = "ph_candidate_value",
  "Lactate" = "lactate_candidate_value", "Hemoglobin" = "hb_candidate_value",
  "Concurrent vasoactive medication use" = "concurrent_vasoactive_medication_use_state",
  "Current invasive mechanical ventilation" = "current_invasive_mechanical_ventilation_state")
if (any(!summary$feature_group %in% names(feature_sources)) ||
    any(!unname(feature_sources[summary$feature_group]) %in% names(values))) {
  stop("RF_SHAP_FEATURE_VALUE_MAPPING_MISSING", call. = FALSE)
}
feature_order <- summary$feature_group[order(summary$mean_abs_shap, decreasing = TRUE)]
categorical <- c("Sex", "Concurrent vasoactive medication use",
                 "Current invasive mechanical ventilation")
long <- rbindlist(lapply(summary$feature_group, function(feature_value) {
  j <- match(feature_value, colnames(grouped))
  source_name <- unname(feature_sources[[feature_value]])
  raw <- values[[source_name]]
  if (identical(source_name, "age_interval")) {
    raw <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", as.character(raw))))
  }
  score <- rep(NA_real_, length(raw))
  if (!feature_value %in% categorical) {
    numeric <- suppressWarnings(as.numeric(raw))
    finite <- is.finite(numeric)
    if (sum(finite) > 1L) score[finite] <-
      (rank(numeric[finite], ties.method = "average") - 1) / (sum(finite) - 1)
  }
  data.table(feature_group = feature_value, case_index = seq_len(nrow(grouped)),
    shap_value = grouped[, j], color_score = score)
}))
long[, feature_group := factor(feature_group, levels = rev(feature_order))]
make_dep <- function(group, field) {
  j <- match(group, colnames(grouped))
  x <- suppressWarnings(as.numeric(values[[field]]))
  z <- data.table(feature_value = x, shap_value = grouped[, j])
  z[is.finite(feature_value) & is.finite(shap_value)]
}
dependence <- list(PaO2 = make_dep("PaO2", "pao2_candidate_value"),
                   pH = make_dep("pH", "ph_candidate_value"),
                   PaCO2 = make_dep("PaCO2", "paco2_candidate_value"))
if (any(vapply(dependence, nrow, integer(1L)) < 10L)) {
  stop("RF_SHAP_DEPENDENCE_OBSERVED_N_INSUFFICIENT", call. = FALSE)
}
theme_rf <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 9.5, hjust = 0),
        plot.margin = margin(5, 7, 5, 5))
plot_a <- ggplot(long, aes(x = shap_value, y = feature_group, color = color_score)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0, height = 0.18, seed = 40507),
             size = 0.7, alpha = 0.60) +
  scale_color_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
    limits = c(0, 1), na.value = "grey65", breaks = c(0, 1), labels = c("Low", "High"),
    name = "Feature value") +
  labs(x = "SHAP value for predicted Bias, percentage points", y = NULL,
       title = "SHAP summary") + theme_rf +
  theme(legend.position = "bottom", legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5))
plot_dep <- function(data, title, x_label) {
  ggplot(data, aes(x = feature_value, y = shap_value)) +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_point(color = "black", size = 0.75, alpha = 0.35) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                color = "black", linewidth = 0.7, span = 0.75) +
    labs(x = x_label, y = "SHAP value for predicted Bias, percentage points",
         title = title) + theme_rf
}
plot_b <- plot_dep(dependence$PaO2, "PaO2", "PaO2, mm Hg")
plot_c <- plot_dep(dependence$pH, "pH", "pH")
plot_d <- plot_dep(dependence$PaCO2, "PaCO2", "PaCO2, mm Hg")
figure <- ((plot_a | plot_b) / (patchwork::free(plot_c, type = "label", side = "l") | plot_d)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11),
        plot.tag.position = c(0.01, 0.99))
ggsave(pdf_path, figure, width = 11.69, height = 8.27, units = "in",
       device = cairo_pdf, bg = "white")
receipt <- data.table(figure_id = paste0("Supplementary Figure ", figure_number),
  cohort = cohort, source_mode = status, shap_patient_n = nrow(sample),
  shap_feature_n = ncol(grouped), shap_point_n = nrow(long),
  pao2_observed_n = nrow(dependence$PaO2), ph_observed_n = nrow(dependence$pH),
  paco2_observed_n = nrow(dependence$PaCO2))
fwrite(receipt, receipt_path, sep = "\t")
cat("RF_SHAP_FIGURE_CANDIDATE_PASS figure=", figure_number,
    " patients=", nrow(sample), "\n", sep = "")
