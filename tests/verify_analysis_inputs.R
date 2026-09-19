args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: ANALYSIS_OUTPUT_DIR NEW_QA_DIR", call. = FALSE)
input <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (file.exists(qa) || dir.exists(qa)) stop("QA_OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
suppressPackageStartupMessages(library(digest))
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
checks <- list()
check <- function(name, value) {
  pass <- isTRUE(value)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = pass)
}
key <- function(x, modality = NULL) {
  fields <- c("patient_key_internal", "stay_key_internal")
  if (!is.null(modality)) fields <- c(fields, paste0(modality, "_stable_event_ordinal"))
  do.call(paste, c(x[fields], sep = "\r"))
}
check("builder_status", identical(readLines(file.path(input, "status.txt"), warn = FALSE),
                                  "ANALYSIS_INPUT_BUILD_PASS"))
receipt <- read.delim(file.path(input, "receipt.tsv"), check.names = FALSE)
expected <- c("counts.tsv", paste0("matched_pairs_", c(60, 5), ".rds"),
              paste0("day1_7_pairs_", c(60, 5), ".rds"),
              paste0("analytic_pairs_", c(60, 5), ".rds"),
              paste0("lmm_input_", c(60, 5), ".rds"), "rf_candidate_60.rds")
check("receipt_complete", setequal(receipt$file, expected) && !anyDuplicated(receipt$file))
check("receipt_hashes", all(vapply(seq_len(nrow(receipt)), function(i) {
  path <- file.path(input, receipt$file[[i]])
  file.exists(path) && identical(tolower(receipt$sha256[[i]]),
                                 digest::digest(path, file = TRUE, algo = "sha256"))
}, logical(1))))
counts <- read.delim(file.path(input, "counts.tsv"), check.names = FALSE)
check("two_independent_windows", identical(sort(as.integer(counts$window_minutes)), c(5L, 60L)))
for (window in c(60L, 5L)) {
  pair <- readRDS(file.path(input, paste0("matched_pairs_", window, ".rds")))
  q21 <- readRDS(file.path(input, paste0("day1_7_pairs_", window, ".rds")))
  analytic <- readRDS(file.path(input, paste0("analytic_pairs_", window, ".rds")))
  lmm <- readRDS(file.path(input, paste0("lmm_input_", window, ".rds")))
  row <- counts[counts$window_minutes == window, , drop = FALSE]
  check(paste0("row_counts_", window), nrow(row) == 1L && nrow(pair) == nrow(analytic) &&
          nrow(pair) == nrow(lmm) && nrow(pair) == row$pair_n[[1L]])
  check(paste0("day1_7_count_", window), nrow(q21) == row$day1_7_pair_n[[1L]])
  check(paste0("day1_7_scope_", window), !nrow(q21) ||
          all(q21$spo2_relative_icu_minutes >= 0 & q21$sao2_relative_icu_minutes >= 0 &
              q21$spo2_relative_icu_minutes <= 10080 & q21$sao2_relative_icu_minutes <= 10080))
  check(paste0("day1_7_one_to_one_", window), !anyDuplicated(key(q21, "spo2")) &&
          !anyDuplicated(key(q21, "sao2")))
  check(paste0("day1_7_labels_", window), !nrow(q21) ||
          all(q21$relative_icu_day == pmin(7L, floor(q21$sao2_relative_icu_minutes / 1440) + 1L) &
              q21$sao2_lt88 == (q21$sao2_saturation_percent < 88) &
              q21$spo2_display_stratum == ifelse(q21$spo2_saturation_percent < 88, "70_to_lt88",
                 ifelse(q21$spo2_saturation_percent < 92, "88_to_lt92",
                        ifelse(q21$spo2_saturation_percent <= 96, "92_to_96", "gt96_to_100")))))
  check(paste0("one_to_one_", window), !anyDuplicated(key(pair, "spo2")) &&
          !anyDuplicated(key(pair, "sao2")))
  check(paste0("lag_", window), all(is.finite(pair$absolute_lag_minutes)) &&
          all(abs(pair$spo2_relative_icu_minutes - pair$sao2_relative_icu_minutes) ==
              pair$absolute_lag_minutes) && all(pair$absolute_lag_minutes <= window))
  check(paste0("identity_", window), identical(as.character(pair$patient_key_internal),
            as.character(analytic$patient_key_internal)) &&
          identical(as.character(pair$stay_key_internal), as.character(analytic$stay_key_internal)) &&
          identical(as.character(analytic$pair_key_internal), as.character(lmm$pair_key_internal)))
  check(paste0("bias_", window), isTRUE(all.equal(analytic$bias_spo2_minus_sao2,
          pair$spo2_saturation_percent - pair$sao2_saturation_percent,
          check.attributes = FALSE)))
  check(paste0("paired_mean_", window), isTRUE(all.equal(
          lmm$paired_mean_saturation_c90,
          (pair$spo2_saturation_percent + pair$sao2_saturation_percent) / 2 - 90,
          check.attributes = FALSE)))
  bounds <- list(pao2 = c(0, 1000), paco2 = c(0, 1000),
                 ph = c(6.3, 8.0), lactate = c(0, Inf), hb = c(0, 30))
  valid_values <- vapply(names(bounds), function(metric) {
    value <- analytic[[paste0(metric, "_candidate_value")]]
    observed <- value[!is.na(value)]
    lower <- bounds[[metric]][[1L]]
    upper <- bounds[[metric]][[2L]]
    if (metric == "ph") all(observed >= lower & observed <= upper) else
      all(observed > lower & observed <= upper)
  }, logical(1))
  check(paste0("clinical_validity_bounds_", window), all(valid_values))
  check(paste0("clinical_validity_state_", window), all(vapply(names(bounds), function(metric) {
    value <- analytic[[paste0(metric, "_candidate_value")]]
    state <- analytic[[paste0(metric, "_value_state")]]
    all(is.na(value) | state == "eligible_pending_candidate_build")
  }, logical(1))))
  check(paste0("anchor_", window), isTRUE(all.equal(lmm$sao2_anchor_icu_hours,
          pair$sao2_relative_icu_minutes / 60, check.attributes = FALSE)))
  check(paste0("count_units_", window),
        length(unique(key(pair))) == row$stay_n[[1L]] &&
          length(unique(as.character(pair$patient_key_internal))) == row$patient_n[[1L]])
}
rf <- readRDS(file.path(input, "rf_candidate_60.rds"))
analytic <- readRDS(file.path(input, "analytic_pairs_60.rds"))
check("rf_identity", identical(as.character(rf$pair_uid), as.character(analytic$pair_key_internal)))
check("rf_outcome", isTRUE(all.equal(rf$bias_spo2_minus_sao2,
       analytic$bias_spo2_minus_sao2, check.attributes = FALSE)))
check("rf_paired_mean", isTRUE(all.equal(rf$paired_mean_saturation_c90,
      (rf$spo2_percent + rf$sao2_percent) / 2 - 90, check.attributes = FALSE)))
check("rf_age_representation", if (counts$cohort[[1L]] %in% c("MIMIC", "eICU"))
      all(!is.na(rf$age_topcoded_indicator)) else TRUE)
table <- do.call(rbind, checks)
utils::write.table(table, file.path(qa, "checks.tsv"), sep = "\t", row.names = FALSE,
                   quote = FALSE)
status <- if (all(table$pass)) "ANALYSIS_INPUT_INDEPENDENT_QA_PASS" else
  "ANALYSIS_INPUT_INDEPENDENT_QA_FAIL"
writeLines(status, file.path(qa, "status.txt"))
cat(status, "checks=", nrow(table), "\n")
if (!all(table$pass)) quit(save = "no", status = 1L)
