#!/usr/bin/env Rscript
# Supplementary Figure 4: mean Bias points and 95% LoA across SpO2 strata.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: MANIFEST.tsv NEW_OUTPUT.pdf NEW_RECEIPT.tsv", call. = FALSE)
manifest <- fread(args[[1L]])
pdf_path <- args[[2L]]
receipt_path <- args[[3L]]
if (any(file.exists(c(pdf_path, receipt_path)))) stop("OUTPUT_EXISTS_REFUSING_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(c(pdf_path, receipt_path))))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
if (!all(c("cohort", "agreement_tsv") %in% names(manifest))) stop("MANIFEST_COLUMNS_INVALID", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
labels <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
if (nrow(manifest) != 5L || !setequal(manifest$cohort, cohorts) || anyDuplicated(manifest$cohort)) {
  stop("FIVE_COHORT_SET_INVALID", call. = FALSE)
}
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
required <- c("cohort", "window_minutes", "spo2_stratum", "pair_n", "patient_n",
              "mean_bias", "lower_loa", "upper_loa")
data <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  source <- file.path(dirname(manifest$agreement_tsv[[i]]), "agreement_by_spo2_stratum.tsv")
  if (!file.exists(source)) stop("STRATUM_AGREEMENT_INPUT_MISSING", call. = FALSE)
  x <- fread(source)
  if (!all(required %in% names(x))) stop("STRATUM_AGREEMENT_COLUMNS_MISSING", call. = FALSE)
  if (nrow(x) != 4L || any(x$cohort != manifest$cohort[[i]]) ||
      any(x$window_minutes != 60L) || !setequal(x$spo2_stratum, strata)) {
    stop(paste("STRATUM_AGREEMENT_GRID_INVALID", manifest$cohort[[i]]), call. = FALSE)
  }
  x
}), use.names = TRUE, fill = TRUE)
if (nrow(data) != 20L || anyDuplicated(data[, paste(cohort, spo2_stratum)]) ||
    anyNA(data[, ..required]) ||
    any(data$pair_n < 1 | data$patient_n < 1 | data$patient_n > data$pair_n |
        !is.finite(data$lower_loa) | !is.finite(data$mean_bias) |
        !is.finite(data$upper_loa) | data$lower_loa >= data$mean_bias |
        data$mean_bias >= data$upper_loa)) {
  stop("STRATUM_AGREEMENT_VALUES_INVALID", call. = FALSE)
}
data[, `:=`(cohort_label = factor(unname(labels[cohort]), levels = unname(labels[cohorts])),
             stratum_label = factor(spo2_stratum, levels = strata,
                                    labels = c("70%\u2013<88%", "88%\u2013<92%",
                                               "92%\u2013\u226496%", ">96%\u2013100%")))]
x_limit <- max(30, ceiling(max(abs(c(data$lower_loa,data$upper_loa)))/5)*5)
if (!is.finite(x_limit) || x_limit > 1000)
  stop("STRATUM_AGREEMENT_AXIS_RANGE_INVALID", call. = FALSE)
p <- ggplot(data, aes(y = stratum_label, x = mean_bias)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.25) +
  geom_segment(aes(x = lower_loa, xend = upper_loa, yend = stratum_label), linewidth = 0.5) +
  geom_errorbar(aes(xmin = lower_loa, xmax = upper_loa), orientation = "y",
                width = 0.17, linewidth = 0.35) +
  geom_point(size = 1.9) +
  facet_wrap(~ cohort_label, ncol = 3, drop = FALSE) +
  scale_x_continuous(limits = c(-x_limit, x_limit),
                     breaks = pretty(c(-x_limit, x_limit), n = 7)) +
  scale_y_discrete(limits = levels(data$stratum_label)) +
  labs(x = "Mean Bias (point) and 95% LoA (line), percentage points",
       y = expression(SpO[2]~stratum)) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
                                  strip.background = element_blank(),
                                  strip.text = element_text(face = "bold"))
ggsave(pdf_path, p, width = 12, height = 7.8, device = cairo_pdf, bg = "white")
receipt <- data[, .(stratum_n = .N, pair_n_total = sum(pair_n)), by = .(cohort, window_minutes)]
receipt[, x_axis_abs_limit := x_limit]
fwrite(receipt, receipt_path, sep = "\t")
cat("STRATUM_AGREEMENT_CANDIDATE_RENDER_PASS cells=20\n")
