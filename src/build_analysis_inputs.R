args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: COHORT STANDARDIZED_INPUT_DIR NEW_OUTPUT_DIR", call. = FALSE)
cohort <- args[[1L]]
input <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (!cohort %in% c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang"))
  stop("COHORT_INVALID", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages({library(data.table); library(digest)})
for (name in c("r04_event_core.R", "r05_pairing_core.R", "r211b_core.R"))
  source(file.path(root, "code", "cores", name), local = FALSE)
source(file.path(root, "code", "cores", "stage404_core.R"), local = FALSE)
source(file.path(root, "code", "formal", "stage405_multilevel_core.R"), local = FALSE)
source(file.path(root, "code", "formal", "stage4_variable_validity_core.R"), local = FALSE)
for (name in c("standardized_event_input_v1.R", "standardized_clinical_input_v1.R",
               "model_input_adapter_v1.R"))
  source(file.path(root, "code", name), local = FALSE)
required <- c("stays.tsv", "spo2_events.tsv", "sao2_events.tsv", "stay_clinical.tsv",
              "pair_clinical_60.tsv", "pair_clinical_5.tsv")
if (!all(file.exists(file.path(input, required))))
  stop(paste("STANDARDIZED_INPUT_MISSING:", paste(required[!file.exists(file.path(input, required))], collapse = ", ")), call. = FALSE)
loaded <- si_load_input(input, cohort)
stays <- ci1_read_stays(file.path(input, "stay_clinical.tsv"), loaded$stays)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output)) stop("OUTPUT_CREATE_FAILED", call. = FALSE)
writeLines("ANALYSIS_INPUT_BUILD_RUNNING", file.path(output, "status.txt"))
tryCatch({
  rows <- list()
  for (window in c(60L, 5L)) {
    matched <- r05_pair_cohort(loaded$spo2, loaded$sao2, cohort, window)$final_pairs
    if (!nrow(matched)) stop(paste0("NO_ELIGIBLE_PAIRS_", window), call. = FALSE)
    q21_spo2 <- loaded$spo2[loaded$spo2$event_time_min <= 10080, , drop = FALSE]
    q21_sao2 <- loaded$sao2[loaded$sao2$event_time_min <= 10080, , drop = FALSE]
    q21 <- if (nrow(q21_spo2) && nrow(q21_sao2))
      r05_pair_cohort(q21_spo2, q21_sao2, cohort, window)$final_pairs else matched[0, ]
    if (nrow(q21)) {
      q21$relative_icu_day <- pmin(7L, floor(q21$sao2_relative_icu_minutes / 1440) + 1L)
      q21$spo2_display_stratum <- s404_source_stratum_code(q21$spo2_saturation_percent)
      q21$sao2_lt88 <- q21$sao2_saturation_percent < 88
      s404_validate_q21_pairs(data.table::as.data.table(q21), cohort,
                              if (window == 60L) "M60" else "M5")
    }
    saveRDS(q21, file.path(output, paste0("day1_7_pairs_", window, ".rds")))
    concurrent <- ci1_read_pairs(file.path(input, paste0("pair_clinical_", window, ".tsv")),
                                 matched, stays, cohort, window)
    analytic <- ci1_make_analytic_pairs(matched, stays, concurrent, cohort, window)
    metric_map <- c(pao2 = "pao2_mmhg", paco2 = "paco2_mmhg", ph = "ph",
                    lactate = "lactate_mmol_l", hb = "hb_g_dl")
    for (metric in names(metric_map)) {
      field <- paste0(metric, "_candidate_value")
      state <- paste0(metric, "_value_state")
      checked <- s4v_validate(analytic[[field]], metric_map[[metric]], analytic[[state]])
      analytic[[field]] <- checked$value
      analytic[[state]] <- checked$value_state
      analytic[[paste0(metric, "_diagnostic_review")]] <- checked$diagnostic_review
    }
    lmm <- mi1_lmm_frame(analytic, matched, cohort, window)
    if (nrow(analytic) != nrow(matched) || nrow(lmm) != nrow(matched))
      stop("ANALYSIS_INPUT_ROW_COUNT_MISMATCH", call. = FALSE)
    saveRDS(matched, file.path(output, paste0("matched_pairs_", window, ".rds")))
    saveRDS(analytic, file.path(output, paste0("analytic_pairs_", window, ".rds")))
    saveRDS(lmm, file.path(output, paste0("lmm_input_", window, ".rds")))
    if (window == 60L) {
      rf <- mi1_rf_candidate(analytic, cohort)
      rf$paired_mean_saturation_c90 <-
        (rf$spo2_percent + rf$sao2_percent) / 2 - 90
      rf <- as.data.frame(s405m_derive_age_fields(data.table::as.data.table(rf), cohort))
      if (nrow(rf) != nrow(matched)) stop("RF_INPUT_ROW_COUNT_MISMATCH", call. = FALSE)
      saveRDS(rf, file.path(output, "rf_candidate_60.rds"))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = cohort, window_minutes = window, pair_n = nrow(matched),
      stay_n = length(unique(as.character(analytic$stay_key_internal))),
      patient_n = length(unique(as.character(analytic$patient_key_internal))),
      day1_7_pair_n = nrow(q21))
  }
  utils::write.table(do.call(rbind, rows), file.path(output, "counts.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  files <- list.files(output, pattern = "\\.rds$|^counts\\.tsv$", full.names = TRUE)
  receipt <- data.frame(file = basename(files),
                        sha256 = vapply(files, digest::digest, character(1),
                                        file = TRUE, algo = "sha256"),
                        stringsAsFactors = FALSE)
  utils::write.table(receipt, file.path(output, "receipt.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines("ANALYSIS_INPUT_BUILD_PASS", file.path(output, "status.txt"))
  cat("ANALYSIS_INPUT_BUILD_PASS\n")
}, error = function(e) {
  writeLines(paste("ANALYSIS_INPUT_BUILD_FAIL", conditionMessage(e)),
             file.path(output, "status.txt"))
  stop(e)
})
