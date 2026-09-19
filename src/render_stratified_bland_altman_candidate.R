#!/usr/bin/env Rscript
# Supplementary Figure 3: 60-minute Bland-Altman panels by displayed SpO2 stratum.
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
script_flag <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_flag) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)))
source(file.path(root, "code", "cores", "agreement_core.R"), local = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: MANIFEST60.tsv NEW.pdf NEW_RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
pdf_path <- args[[2L]]
receipt_path <- args[[3L]]
if (file.exists(pdf_path) || file.exists(receipt_path)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (!dir.exists(dirname(pdf_path)) || !dir.exists(dirname(receipt_path))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
strata_labels <- c("70%\u2013<88%", "88%\u2013<92%", "92%\u2013\u226496%", ">96%\u2013100%")
if (!identical(names(m), c("cohort", "analysis_rds", "agreement_tsv")) ||
    nrow(m) != 5L || !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort)) {
  stop("STRATIFIED_BA_MANIFEST_INVALID", call. = FALSE)
}
points <- vector("list", 5L)
smooth <- vector("list", 5L)
lines <- vector("list", 5L)
receipt <- vector("list", 5L)
for (i in seq_len(nrow(m))) {
  cohort_value <- m$cohort[[i]]
  stratum_path <- file.path(dirname(m$agreement_tsv[[i]]), "agreement_by_spo2_stratum.tsv")
  if (!file.exists(m$analysis_rds[[i]]) || !file.exists(stratum_path)) {
    stop("STRATIFIED_BA_INPUT_MISSING", call. = FALSE)
  }
  x <- as.data.table(readRDS(m$analysis_rds[[i]]))
  s <- fread(stratum_path)
  required <- c("cohort", "spo2_saturation_percent", "sao2_saturation_percent",
                "bias_spo2_minus_sao2", "pair_window_minutes")
  if (!all(required %in% names(x)) || nrow(x) < 1L ||
      !all(x$cohort == cohort_value) || !all(x$pair_window_minutes == 60L) ||
      nrow(s) != 4L || !setequal(s$spo2_stratum, strata) ||
      !all(s$cohort == cohort_value) || !all(s$window_minutes == 60L)) {
    stop("STRATIFIED_BA_SOURCE_CONTRACT_INVALID", call. = FALSE)
  }
  x[, paired_mean := (spo2_saturation_percent + sao2_saturation_percent) / 2]
  x[, bias := spo2_saturation_percent - sao2_saturation_percent]
  if (any(!is.finite(x$bias)) || any(abs(x$bias - x$bias_spo2_minus_sao2) > 1e-8)) {
    stop("STRATIFIED_BA_BIAS_ARITHMETIC_INVALID", call. = FALSE)
  }
  x[, spo2_stratum := fifelse(spo2_saturation_percent < 88, strata[[1L]],
      fifelse(spo2_saturation_percent < 92, strata[[2L]],
        fifelse(spo2_saturation_percent <= 96, strata[[3L]], strata[[4L]])))]
  for (stratum_value in strata) {
    y <- x[spo2_stratum == stratum_value]
    estimate <- s[spo2_stratum == stratum_value]
    if (nrow(estimate) != 1L || estimate$pair_n[[1L]] != nrow(y) ||
        !is.finite(estimate$mean_bias[[1L]]) ||
        estimate$lower_loa[[1L]] >= estimate$mean_bias[[1L]] ||
        estimate$upper_loa[[1L]] <= estimate$mean_bias[[1L]]) {
      stop(paste("STRATIFIED_BA_AGGREGATE_MISMATCH", cohort_value, stratum_value), call. = FALSE)
    }
  }
  points[[i]] <- x[, .(cohort = cohort_value, spo2_stratum, paired_mean, bias)]
  smooth[[i]] <- rbindlist(lapply(strata, function(stratum_value) {
    y <- points[[i]][spo2_stratum == stratum_value]
    curve <- as.data.table(descriptive_smooth(copy(y)))
    curve[, `:=`(cohort = cohort_value, spo2_stratum = stratum_value)]
    curve
  }))
  lines[[i]] <- rbindlist(lapply(strata, function(stratum_value) {
    estimate <- s[spo2_stratum == stratum_value]
    data.table(cohort = cohort_value, spo2_stratum = stratum_value,
      line_type = c("Mean Bias", "Lower 95% LoA", "Upper 95% LoA"),
      value = c(estimate$mean_bias, estimate$lower_loa, estimate$upper_loa))
  }))
  receipt[[i]] <- s[, .(cohort, window_minutes, spo2_stratum, pair_n, mean_bias,
                       lower_loa, upper_loa)]
}
point_data <- rbindlist(points)
smooth_data <- rbindlist(smooth)
line_data <- rbindlist(lines)
receipts <- rbindlist(receipt)
for (dataset in c("point_data", "smooth_data", "line_data")) {
  value <- get(dataset)
  value[, cohort := factor(cohort, levels = cohorts, labels = pretty)]
  value[, spo2_stratum := factor(spo2_stratum, levels = strata, labels = strata_labels)]
  assign(dataset, value)
}
line_data[, line_type := factor(line_type,
  levels = c("Mean Bias", "Lower 95% LoA", "Upper 95% LoA"))]
plot_page <- function(selected) {
  ggplot(point_data[cohort %in% selected], aes(x = paired_mean, y = bias)) +
    geom_point(color = "black", alpha = 0.12, size = 0.12) +
    geom_hline(data = line_data[cohort %in% selected],
      aes(yintercept = value, linetype = line_type), color = "black",
      linewidth = 0.32, inherit.aes = FALSE) +
    geom_line(data = smooth_data[cohort %in% selected],
      aes(x = paired_mean, y = smooth_bias), color = "grey50",
      linewidth = 0.52, inherit.aes = FALSE) +
    facet_grid(cohort ~ spo2_stratum, drop = TRUE, switch = "y") +
    scale_linetype_manual(values = c("solid", "dashed", "dashed"), guide = "none") +
    scale_x_continuous(limits = c(70, 100), breaks = c(70, 80, 90, 100)) +
    scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, 10)) +
    labs(x = "Paired mean saturation, %", y = "Bias (SpO2 - SaO2), percentage points") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.spacing = grid::unit(0.1, "cm"),
          strip.background = element_blank(), strip.placement = "outside",
          strip.text = element_text(face = "bold"))
}
grDevices::cairo_pdf(pdf_path, width = 11.69, height = 8.27, family = "sans", bg = "white")
print(plot_page(unname(pretty[cohorts[1:3]])))
print(plot_page(unname(pretty[cohorts[4:5]])))
dev.off()
fwrite(receipts, receipt_path, sep = "\t")
cat("STRATIFIED_BA_CANDIDATE_RENDER_PASS facets=20 points=", nrow(point_data), "\n", sep = "")
