#!/usr/bin/env Rscript
# Invented seven-day aggregate fixture for display-layout tests only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: NEW_DIRECTORY WINDOW(60|5)", call. = FALSE)
out <- args[[1L]]
window <- suppressWarnings(as.integer(args[[2L]]))
if (is.na(window) || !(window %in% c(5L, 60L))) stop("WINDOW_INVALID", call. = FALSE)
if (file.exists(out)) stop("FIXTURE_OUTPUT_EXISTS", call. = FALSE)
dir.create(out, recursive = TRUE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
manifest <- data.frame(cohort = cohorts, day_tsv = character(length(cohorts)))
for (i in seq_along(cohorts)) {
  rows <- expand.grid(relative_icu_day = 1:7, spo2_stratum = strata,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  rows <- rows[order(rows$relative_icu_day, match(rows$spo2_stratum, strata)), ]
  base <- c(0.44, 0.17, 0.035, 0.015)[match(rows$spo2_stratum, strata)]
  variation <- (i - 3) * 0.006 + sin(rows$relative_icu_day * 0.9) * 0.012
  nominal <- pmax(0.003, pmin(0.9, base + variation))
  rows$denominator_n <- 500L
  rows$patient_n <- 150L
  sparse <- cohorts[[i]] == "Lianyungang" &
    ((rows$relative_icu_day == 5L & rows$spo2_stratum == strata[[3L]]) |
     (rows$relative_icu_day == 6L & rows$spo2_stratum == strata[[4L]]))
  rows$denominator_n[sparse] <- 12L
  rows$patient_n[sparse] <- 6L
  rows$numerator_n <- as.integer(round(rows$denominator_n * nominal))
  rows$conditional_proportion <- rows$numerator_n / rows$denominator_n
  half_width <- ifelse(rows$spo2_stratum %in% strata[1:2], 0.035, 0.008)
  rows$ci_lower <- pmax(0, rows$conditional_proportion - half_width)
  rows$ci_upper <- pmin(1, rows$conditional_proportion + half_width)
  rows$cell_status <- ifelse(sparse, "SPARSE", "DISPLAYABLE")
  rows$conditional_proportion[sparse] <- NA_real_
  rows$ci_lower[sparse] <- NA_real_
  rows$ci_upper[sparse] <- NA_real_
  rows$cohort <- cohorts[[i]]
  rows$window_minutes <- window
  path <- normalizePath(file.path(out, paste0(cohorts[[i]], "_day.tsv")), mustWork = FALSE)
  write.table(rows, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  manifest$day_tsv[[i]] <- path
}
write.table(manifest, file.path(out, "manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("SEVEN_DAY_INVENTED_FIXTURE_PASS rows=140 sparse=2 window=", window, "\n", sep = "")
