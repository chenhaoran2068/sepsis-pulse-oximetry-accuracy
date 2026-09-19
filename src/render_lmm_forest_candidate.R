#!/usr/bin/env Rscript
# Five-cohort LMM coefficient forest from the formal pooled, unscaled outputs.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: MANIFEST.tsv WINDOW(60|5) NEW_OUTPUT.pdf NEW_CHART_DATA.tsv NEW_RECEIPT.tsv", call. = FALSE)
}
m <- fread(args[[1L]])
window <- suppressWarnings(as.integer(args[[2L]]))
pdf_path <- args[[3L]]
data_path <- args[[4L]]
receipt_path <- args[[5L]]
if (is.na(window) || !window %in% c(60L, 5L)) stop("WINDOW_INVALID", call. = FALSE)
if (any(file.exists(c(pdf_path, data_path, receipt_path)))) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(c(pdf_path, data_path, receipt_path))))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
if (!identical(names(m), c("cohort", "lmm_tsv"))) stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
labels <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
if (nrow(m) != 5L || !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort)) {
  stop("FIVE_COHORT_SET_INVALID", call. = FALSE)
}
if (any(!file.exists(m$lmm_tsv))) stop("MANIFEST_INPUT_MISSING", call. = FALSE)
variable_map <- data.table(
  variable_id = c("age", "sex", "cci", "sofa", "pao2", "paco2", "ph", "lactate", "hb", "paired_mean"),
  term = c("age_value_per_10y", "sex_classfemale", "cci_value", "sofa_total",
           "pao2_candidate_value", "paco2_candidate_value", "ph_candidate_value",
           "lactate_candidate_value", "hb_candidate_value", "paired_mean_saturation_c90"),
  display = c("Age, per 10 y", "Female vs male", "Charlson Comorbidity Index, per point",
              "SOFA, per point", "PaO2, per 10 mm Hg", "PaCO2, per 10 mm Hg", "pH, per 0.10",
              "Lactate, per 1 mmol/L", "Hemoglobin, per 1 g/dL",
              "Paired mean saturation, per 1 percentage point"),
  display_order = 1:10
)
required <- c("term", "estimate", "ci_lower", "ci_upper", "reporting_multiplier", "imputations")
all <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
  x <- fread(m$lmm_tsv[[i]])
  if (!all(required %in% names(x))) stop("LMM_COLUMNS_MISSING", call. = FALSE)
  if (anyDuplicated(x$term)) stop("LMM_DUPLICATE_TERM", call. = FALSE)
  x[, cohort := m$cohort[[i]]]
  x
}), use.names = TRUE, fill = TRUE)
if (anyNA(all[, ..required]) || any(all$reporting_multiplier <= 0) ||
    any(all$imputations != 30L)) stop("LMM_INPUT_VALUE_INVALID", call. = FALSE)
direct <- merge(all, variable_map, by = "term", all = FALSE)
if (any(direct$ci_lower > direct$estimate | direct$estimate > direct$ci_upper)) {
  stop("LMM_CI_INVALID", call. = FALSE)
}
direct[, `:=`(estimate_scaled = estimate * reporting_multiplier,
             lower_scaled = ci_lower * reporting_multiplier,
             upper_scaled = ci_upper * reporting_multiplier)]
grid <- CJ(cohort = cohorts, variable_id = variable_map$variable_id)
grid <- merge(grid, variable_map, by = "variable_id", all.x = TRUE, sort = FALSE)
grid <- merge(grid, direct[, .(cohort, variable_id, estimate_scaled, lower_scaled,
                              upper_scaled, reporting_multiplier)],
              by = c("cohort", "variable_id"), all.x = TRUE, sort = FALSE)
if (nrow(grid) != 50L || anyDuplicated(grid[, paste(cohort, variable_id)])) {
  stop("LMM_DISPLAY_GRID_INVALID", call. = FALSE)
}
grid[, `:=`(status = fifelse(is.finite(estimate_scaled), "estimated",
                           fifelse(variable_id == "age" & cohort %chin% c("Amsterdam", "SICDB"),
                                   "categorical_age", "not_in_model")),
             cohort_label = factor(unname(labels[cohort]), levels = unname(labels[cohorts])),
             y_position = 11L - display_order)]
if (grid[variable_id == "paired_mean", any(status != "estimated")]) {
  stop("PAIRED_MEAN_MISSING", call. = FALSE)
}
signed <- function(z) {
  z[abs(z) < 0.005] <- 0
  sprintf("%+.2f", z)
}
grid[, estimate_label := fifelse(status == "estimated",
                                paste0(signed(estimate_scaled), " (", signed(lower_scaled),
                                       " to ", signed(upper_scaled), ")"),
                                fifelse(status == "categorical_age", "Categorical age (Table S20)", "-"))]
grid <- grid[order(match(cohort, cohorts), display_order)]
est <- grid[status == "estimated"]
effect_limit <- ceiling(max(abs(c(est$lower_scaled, est$upper_scaled))) * 10) / 10 + 0.05
if (!is.finite(effect_limit) || effect_limit <= 0) stop("LMM_EFFECT_LIMIT_INVALID", call. = FALSE)
text_start <- effect_limit + 0.08
text_end <- text_start + max(1.8, 3 * effect_limit)
header <- unique(grid[, .(cohort_label)])
header[, `:=`(x = text_start, y = 10.78, label = "Estimate (95% CI)")]
term_labels <- list(
  expression(Age * ", per 10 y")[[1L]], expression("Female vs male")[[1L]],
  expression("Charlson Comorbidity Index, per point")[[1L]], expression(SOFA * ", per point")[[1L]],
  expression(PaO[2] * ", per 10 mm Hg")[[1L]], expression(PaCO[2] * ", per 10 mm Hg")[[1L]],
  expression(pH * ", per 0.10")[[1L]], expression("Lactate, per 1 mmol/L")[[1L]],
  expression("Hemoglobin, per 1 g/dL")[[1L]],
  expression("Paired mean saturation, per 1 percentage point")[[1L]]
)
p <- ggplot(grid, aes(y = y_position)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey35") +
  geom_vline(xintercept = effect_limit + 0.035, linewidth = 0.3, color = "grey70") +
  geom_errorbar(data = est, aes(x = estimate_scaled, xmin = lower_scaled, xmax = upper_scaled),
                orientation = "y", width = 0.17, linewidth = 0.55, color = "black") +
  geom_point(data = est, aes(x = estimate_scaled), shape = 21, fill = "white",
             stroke = 0.5, size = 1.6) +
  geom_text(aes(x = text_start, label = estimate_label), hjust = 0, size = 2.45,
            color = "grey20") +
  geom_text(data = header, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, size = 2.5, fontface = "bold") +
  facet_wrap(~ cohort_label, ncol = 3, drop = FALSE) +
  scale_y_continuous(limits = c(0.45, 11.05), breaks = 10:1, labels = term_labels,
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(-effect_limit, text_end),
                     breaks = pretty(c(-effect_limit, effect_limit), n = 5),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = expression("Difference in " * (SpO[2] - SaO[2]) * " Bias, percentage points"), y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.7),
        axis.line.x = element_blank(), panel.spacing.x = unit(10, "pt"),
        panel.spacing.y = unit(9, "pt"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        plot.margin = margin(14, 12, 8, 8))
ggsave(pdf_path, p, width = 13.2, height = 9.7, device = cairo_pdf, bg = "white")
fwrite(grid[, .(cohort, window_minutes = window, variable_id, term, display_order,
                status, estimate_scaled, lower_scaled, upper_scaled, estimate_label)],
       data_path, sep = "\t")
receipt <- grid[, .(variable_n = .N, estimated_n = sum(status == "estimated"),
                    categorical_age_n = sum(status == "categorical_age"),
                    not_in_model_n = sum(status == "not_in_model")), by = cohort]
receipt[, window_minutes := window]
fwrite(receipt, receipt_path, sep = "\t")
cat("LMM_FOREST_CANDIDATE_RENDER_PASS window=", window, " cells=50\n", sep = "")
