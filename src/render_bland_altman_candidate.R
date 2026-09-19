#!/usr/bin/env Rscript
# Candidate renderer for the five-cohort 60- or 5-minute Bland-Altman display.
# It never reads a protected study path implicitly and refuses existing output.

suppressPackageStartupMessages(library(ggplot2))
script_flag <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_flag) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
candidate_root <- dirname(dirname(script_path))
source(file.path(candidate_root, "code", "cores", "agreement_core.R"), local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript src/render_bland_altman_candidate.R MANIFEST.tsv WINDOW(60|5) NEW_OUTPUT.pdf NEW_RECEIPT.tsv", call. = FALSE)
}
manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
window <- suppressWarnings(as.integer(args[[2L]]))
if (is.na(window) || !(window %in% c(60L, 5L))) stop("WINDOW_INVALID", call. = FALSE)
output_path <- args[[3L]]
receipt_path <- args[[4L]]
if (file.exists(output_path) || file.exists(receipt_path)) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (!dir.exists(dirname(output_path)) || !dir.exists(dirname(receipt_path))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)

manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
required_manifest <- c("cohort", "analysis_rds", "agreement_tsv")
if (!identical(names(manifest), required_manifest)) stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (!identical(sort(manifest$cohort), sort(cohorts))) stop("FIVE_COHORT_SET_INVALID", call. = FALSE)
if (any(!file.exists(manifest$analysis_rds)) || any(!file.exists(manifest$agreement_tsv))) {
  stop("MANIFEST_INPUT_MISSING", call. = FALSE)
}

points <- vector("list", nrow(manifest))
lines <- vector("list", nrow(manifest))
smoothers <- vector("list", nrow(manifest))
receipts <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  cohort <- manifest$cohort[[i]]
  x <- readRDS(manifest$analysis_rds[[i]])
  s <- read.delim(manifest$agreement_tsv[[i]], check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("cohort", "spo2_saturation_percent", "sao2_saturation_percent", "bias_spo2_minus_sao2", "pair_window_minutes")
  if (!all(needed %in% names(x))) stop(paste("PAIR_COLUMNS_MISSING", cohort), call. = FALSE)
  if (!all(c("cohort", "window_minutes", "pair_n", "mean_bias", "lower_loa", "upper_loa") %in% names(s)) || nrow(s) != 1L) {
    stop(paste("AGREEMENT_SCHEMA_INVALID", cohort), call. = FALSE)
  }
  if (nrow(x) < 1L || !all(x$cohort == cohort) || !all(x$pair_window_minutes == window) ||
      s$cohort[[1L]] != cohort || s$window_minutes[[1L]] != window || s$pair_n[[1L]] != nrow(x)) {
    stop(paste("COHORT_WINDOW_OR_COUNT_MISMATCH", cohort), call. = FALSE)
  }
  bias <- as.numeric(x$spo2_saturation_percent) - as.numeric(x$sao2_saturation_percent)
  if (any(!is.finite(bias)) || any(abs(bias - x$bias_spo2_minus_sao2) > 1e-8)) {
    stop(paste("PAIR_BIAS_ARITHMETIC_INVALID", cohort), call. = FALSE)
  }
  if (abs(mean(bias) - s$mean_bias[[1L]]) > 1e-5) stop(paste("MEAN_BIAS_MISMATCH", cohort), call. = FALSE)
  if (any(!is.finite(unlist(s[1L, c("mean_bias", "lower_loa", "upper_loa")])))) {
    stop(paste("AGREEMENT_ESTIMATE_NONFINITE", cohort), call. = FALSE)
  }
  if (!(s$lower_loa[[1L]] < s$mean_bias[[1L]] && s$mean_bias[[1L]] < s$upper_loa[[1L]])) {
    stop(paste("LOA_ORDER_INVALID", cohort), call. = FALSE)
  }
  points[[i]] <- data.frame(cohort = cohort,
                            paired_mean = (x$spo2_saturation_percent + x$sao2_saturation_percent) / 2,
                            bias = bias)
  smoothers[[i]] <- as.data.frame(descriptive_smooth(data.table::as.data.table(points[[i]])))
  if (nrow(smoothers[[i]]) < 2L) stop(paste("LOESS_GRID_EMPTY", cohort), call. = FALSE)
  smoothers[[i]]$cohort <- cohort
  lines[[i]] <- data.frame(cohort = cohort,
                           line_type = c("Mean Bias", "Lower 95% LoA", "Upper 95% LoA"),
                           value = c(s$mean_bias[[1L]], s$lower_loa[[1L]], s$upper_loa[[1L]]))
  receipts[[i]] <- data.frame(cohort = cohort, window_minutes = window,
                              pair_n = nrow(x), recomputed_mean_bias = mean(bias),
                              reported_mean_bias = s$mean_bias[[1L]],
                              lower_loa = s$lower_loa[[1L]], upper_loa = s$upper_loa[[1L]])
}

point_data <- do.call(rbind, points)
line_data <- do.call(rbind, lines)
smooth_data <- do.call(rbind, smoothers)
receipt <- do.call(rbind, receipts)
labels <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
point_data$cohort <- factor(point_data$cohort, levels = cohorts, labels = labels)
line_data$cohort <- factor(line_data$cohort, levels = cohorts, labels = labels)
smooth_data$cohort <- factor(smooth_data$cohort, levels = cohorts, labels = labels)
line_data$line_type <- factor(line_data$line_type,
                              levels = c("Mean Bias", "Lower 95% LoA", "Upper 95% LoA"))
plot <- ggplot(point_data, aes(x = paired_mean, y = bias)) +
  geom_point(color = "black", alpha = 0.13, size = 0.18) +
  geom_hline(data = line_data, aes(yintercept = value, linetype = line_type),
             color = "black", linewidth = 0.38, inherit.aes = FALSE) +
  geom_line(data = smooth_data, aes(x = paired_mean, y = smooth_bias),
            color = "gray45", linewidth = 0.65, inherit.aes = FALSE) +
  facet_wrap(~ cohort, ncol = 3, scales = "fixed") +
  scale_linetype_manual(values = c("solid", "dashed", "dashed"), guide = "none") +
  scale_x_continuous(limits = c(70, 100), breaks = seq(70, 100, 5)) +
  scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, 10)) +
  labs(x = expression("Paired mean saturation, " * "%"),
       y = expression("Paired difference (SpO"[2] * " - SaO"[2] * "), percentage points")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.spacing = grid::unit(0.18, "cm"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(output_path, plot = plot, width = 11.7, height = 7.8, units = "in", device = cairo_pdf)
write.table(receipt, receipt_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("BLAND_ALTMAN_CANDIDATE_RENDER_PASS", window, nrow(point_data), "\n")
