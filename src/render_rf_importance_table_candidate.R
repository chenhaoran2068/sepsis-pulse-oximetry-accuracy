#!/usr/bin/env Rscript
# Technical candidate for Supplementary Table 24, grouped mean absolute SHAP.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: RF_MANIFEST.tsv NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
pretty <- c(MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
            SICDB = "SICdb", Lianyungang = "Lianyungang")
features <- c("Paired mean saturation", "Age", "Sex", "Charlson Comorbidity Index",
  "SOFA", "PaO2", "PaCO2", "pH", "Lactate", "Hemoglobin",
  "Concurrent vasoactive medication use", "Current invasive mechanical ventilation")
if (!identical(names(m), c("cohort", "rf_dir")) || nrow(m) != 5L ||
    !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort)) {
  stop("RF_IMPORTANCE_MANIFEST_INVALID", call. = FALSE)
}
outputs <- args[2:4]
if (any(file.exists(outputs))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(outputs)))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
inputs <- lapply(cohorts, function(cohort_value) {
  dir <- m[cohort == cohort_value]$rf_dir
  summary_path <- file.path(dir, "rf_shap_summary.tsv")
  sample_path <- file.path(dir, "rf_shap_sample.tsv")
  if (!file.exists(summary_path) || !file.exists(sample_path)) stop("RF_IMPORTANCE_SOURCE_MISSING", call. = FALSE)
  x <- fread(summary_path)
  sample <- fread(sample_path)
  if (!all(c("feature_group", "feature_role", "mean_abs_shap") %in% names(x)) ||
      !all(c("pair_uid", "analysis_patient_key") %in% names(sample)) ||
      nrow(x) < 1L || nrow(sample) < 1L || nrow(sample) > 1000L ||
      anyDuplicated(x$feature_group) || anyDuplicated(sample$analysis_patient_key) ||
      any(!x$feature_group %in% features) || anyNA(x) ||
      any(!is.finite(x$mean_abs_shap)) || any(x$mean_abs_shap < 0) ||
      sum(x$feature_group == "Paired mean saturation") != 1L ||
      x[feature_group == "Paired mean saturation"]$feature_role != "structural_oxygenation_background" ||
      any(x[feature_group != "Paired mean saturation"]$feature_role != "clinical_characteristic")) {
    stop("RF_IMPORTANCE_SOURCE_CONTRACT_INVALID", call. = FALSE)
  }
  x[, rank := rank(-mean_abs_shap, ties.method = "min")]
  list(values = x, sample_n = nrow(sample))
})
names(inputs) <- cohorts
rows <- data.table(feature = features,
                   role = c("Structural oxygenation background", rep("Clinical characteristic", length(features) - 1L)))
for (cohort_value in cohorts) {
  x <- inputs[[cohort_value]]$values
  rows[[cohort_value]] <- vapply(features, function(feature_value) {
    z <- x[feature_group == feature_value]
    if (nrow(z) == 0L) "-" else sprintf("%.2f (%d)", z$mean_abs_shap, z$rank)
  }, character(1L))
}
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(outputs[[1L]], width = 11.69, height = 8.27, family = "sans", bg = "white")
grid.newpage()
grid.text("Supplementary Table 24. Grouped Random-Forest Feature Importance",
  x = unit(0.42, "in"), y = unit(7.73, "in"), just = c("left", "centre"),
  gp = gpar(fontsize = 10.5, fontface = "bold"))
xs <- c(0.45, 5.2, 6.45, 7.88, 9.05, 10.27)
headers <- c("Feature", vapply(cohorts, function(z) paste0(pretty[[z]], " (n=", inputs[[z]]$sample_n, ")"),
                             character(1L)))
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.94, 6.94), "in"), gp = gpar(lwd = 1.1))
for (j in seq_along(headers)) grid.text(headers[[j]], x = unit(xs[[j]], "in"),
  y = unit(6.66, "in"), just = c("left", "centre"), gp = gpar(fontsize = 7.4, fontface = "bold"))
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(6.41, 6.41), "in"), gp = gpar(lwd = 0.6))
for (i in seq_len(nrow(rows))) {
  y <- 6.14 - (i - 1L) * 0.35
  if (i == 2L) grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(y + 0.18, y + 0.18), "in"),
                          gp = gpar(lwd = 0.4, col = "grey60"))
  values <- c(rows$feature[[i]], vapply(cohorts, function(z) rows[[z]][[i]], character(1L)))
  for (j in seq_along(values)) grid.text(values[[j]], x = unit(xs[[j]], "in"),
    y = unit(y, "in"), just = c("left", "centre"), gp = gpar(fontsize = 7.4))
}
grid.lines(x = unit(c(0.42, 11.27), "in"), y = unit(c(1.87, 1.87), "in"), gp = gpar(lwd = 1.0))
grid.text("Cells are mean absolute SHAP values in Bias percentage points, followed by within-cohort rank.",
  x = unit(0.42, "in"), y = unit(1.55, "in"), just = c("left", "centre"),
  gp = gpar(fontsize = 7.2))
dev.off()
fwrite(rows, outputs[[2L]], sep = "\t")
fwrite(data.table(table_id = "Supplementary Table 24", cohort_n = 5L, feature_n = length(features),
                  sample_n_total = sum(vapply(inputs, `[[`, integer(1L), "sample_n"))),
  outputs[[3L]], sep = "\t")
cat("RF_IMPORTANCE_TABLE_CANDIDATE_PASS features=", nrow(rows), "\n", sep = "")
