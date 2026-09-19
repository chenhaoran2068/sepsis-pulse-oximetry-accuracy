args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || !args[[1L]] %in% c("agreement", "threshold"))
  stop("Usage: agreement|threshold ANALYSIS_INPUT_DIR OUTPUT_DIR NEW_QA_DIR", call. = FALSE)
kind <- args[[1L]]
input <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[[4L]], winslash = "/", mustWork = FALSE)
if (file.exists(qa) || dir.exists(qa)) stop("QA_OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
dir.create(qa, recursive = TRUE)
checks <- list()
check <- function(name, ok) checks[[length(checks) + 1L]] <<-
  data.frame(check = name, pass = isTRUE(ok))
read_tsv <- function(name) read.delim(file.path(output, name), check.names = FALSE)
receipt <- read_tsv("run_receipt.tsv")
window <- receipt$window_minutes[[1L]]
full_path <- file.path(input, paste0("matched_pairs_", window, ".rds"))
full <- readRDS(full_path)
check("cohort_scope", nrow(receipt) == 1L && all(full$cohort == receipt$cohort[[1L]]))
check("pair_window", all(full$pair_window_minutes == window))
check("full_input_hash", tolower(receipt[[if (kind == "agreement") "input_sha256" else "full_sha256"]][[1L]]) ==
        tolower(digest::digest(file = full_path, algo = "sha256")))
check("nonempty_pairs", nrow(full) > 0L)
if (kind == "agreement") {
  overall <- read_tsv("agreement_overall.tsv")
  balanced <- read_tsv("agreement_patient_balanced.tsv")
  strata <- read_tsv("agreement_by_spo2_stratum.tsv")
  proportional <- read_tsv("agreement_proportional_bias.tsv")
  hetero <- read_tsv("agreement_heteroscedasticity.tsv")
  bias <- full$spo2_saturation_percent - full$sao2_saturation_percent
  check("overall_counts", overall$pair_n[[1L]] == nrow(full) &&
        overall$patient_n[[1L]] == length(unique(full$patient_key_internal)))
  check("arms_recomputed", abs(overall$arms[[1L]] - sqrt(mean(bias^2))) < 1e-10)
  patient_mean <- tapply(bias, full$patient_key_internal, mean)
  patient_second <- tapply(bias^2, full$patient_key_internal, mean)
  check("balanced_patient_counts", nrow(balanced) == 1L &&
        balanced$patient_n[[1L]] == length(patient_mean) &&
        balanced$pair_n[[1L]] == nrow(full))
  check("balanced_mean_recomputed",
        abs(balanced$balanced_mean_bias[[1L]] - mean(patient_mean)) < 1e-10)
  check("balanced_arms_recomputed",
        abs(balanced$balanced_arms[[1L]] - sqrt(mean(patient_second))) < 1e-10)
  check("balanced_loa_recomputed", max(abs(
    c(balanced$balanced_lower_loa[[1L]], balanced$balanced_upper_loa[[1L]]) -
      (mean(patient_mean) + c(-1.96, 1.96) *
         sqrt(mean(patient_second) - mean(patient_mean)^2)))) < 1e-8)
  check("loa_recomputed", max(abs(c(overall$lower_loa[[1L]], overall$upper_loa[[1L]]) -
        (overall$mean_bias[[1L]] + c(-1.96, 1.96) * overall$total_difference_sd[[1L]]))) < 1e-8)
  check("overall_intervals_include_estimates",
        overall$mean_bias_ci_lower <= overall$mean_bias & overall$mean_bias <= overall$mean_bias_ci_upper &
        overall$lower_loa_ci_lower <= overall$lower_loa & overall$lower_loa <= overall$lower_loa_ci_upper &
        overall$upper_loa_ci_lower <= overall$upper_loa & overall$upper_loa <= overall$upper_loa_ci_upper)
  expected_stratum <- ifelse(full$spo2_saturation_percent < 88, "70%-<88%",
    ifelse(full$spo2_saturation_percent < 92, "88%-<92%",
      ifelse(full$spo2_saturation_percent <= 96, "92%-<=96%", ">96%-100%")))
  counts <- table(expected_stratum)
  check("strata_cover_all_pairs", nrow(strata) == 4L &&
        sum(strata$pair_n) == nrow(full) &&
        all(strata$pair_n == as.integer(counts[strata$spo2_stratum])))
  check("stratum_intervals_include_estimates",
        all(strata$mean_bias_ci_lower <= strata$mean_bias & strata$mean_bias <= strata$mean_bias_ci_upper))
  check("diagnostic_output_finite", nrow(proportional) == 1L && nrow(hetero) == 1L &&
        is.finite(proportional$slope[[1L]]) && is.finite(hetero$variance_function_parameter[[1L]]))
  check("nested_structure_declared", receipt$nested_random_effects[[1L]] ==
        (receipt$cohort[[1L]] == "eICU"))
} else {
  day_path <- file.path(input, paste0("day1_7_pairs_", window, ".rds"))
  day_pairs <- readRDS(day_path)
  forward <- read_tsv("threshold_from_low_sao2.tsv")
  reverse <- read_tsv("threshold_from_displayed_spo2.tsv")
  patient <- read_tsv("threshold_patient_summary.tsv")
  concentration <- read_tsv("threshold_event_concentration.tsv")
  cells <- read_tsv("threshold_day1_7_strata.tsv")
  check("day_input_hash", tolower(receipt$day_sha256[[1L]]) ==
        tolower(digest::digest(file = day_path, algo = "sha256")))
  check("full_and_day_counts", receipt$full_pair_n[[1L]] == nrow(full) &&
        receipt$day1_7_pair_n[[1L]] == nrow(day_pairs))
  for (threshold in c(88, 92)) {
    f <- forward[forward$threshold_spo2_ge == threshold, ]
    r <- reverse[reverse$displayed_spo2_threshold_ge == threshold, ]
    p <- patient[patient$displayed_spo2_threshold_ge == threshold, ]
    event <- full$sao2_saturation_percent < 88 & full$spo2_saturation_percent >= threshold
    at_risk <- full$sao2_saturation_percent < 88
    displayed <- full$spo2_saturation_percent >= threshold
    check(paste0("forward_counts_", threshold), nrow(f) == 1L &&
          f$pair_denominator_sao2_lt88_n[[1L]] == sum(at_risk) &&
          f$pair_numerator_n[[1L]] == sum(event))
    check(paste0("reverse_counts_", threshold), nrow(r) == 1L &&
          r$conditional_denominator_spo2_ge_n[[1L]] == sum(displayed) &&
          r$conditional_numerator_sao2_lt88_n[[1L]] == sum(event))
    check(paste0("patient_counts_", threshold), nrow(p) == 1L &&
          p$accepted_patient_n[[1L]] == length(unique(full$patient_key_internal)) &&
          p$affected_patient_n[[1L]] == length(unique(full$patient_key_internal[event])))
    check(paste0("concentration_event_count_", threshold),
          concentration$observed_event_pair_n[concentration$displayed_spo2_threshold_ge == threshold] == sum(event))
  }
  check("threshold_nesting", forward$pair_numerator_n[forward$threshold_spo2_ge == 92] <=
        forward$pair_numerator_n[forward$threshold_spo2_ge == 88])
  check("day_cells_complete", nrow(cells) == 28L &&
        nrow(unique(cells[, c("relative_icu_day", "spo2_stratum")])) == 28L)
  check("day_counts_cover_pairs", sum(cells$denominator_n) == nrow(day_pairs) &&
        sum(cells$numerator_n) == sum(day_pairs$sao2_saturation_percent < 88))
  sparse <- cells$denominator_n < 20L | cells$patient_n < 10L
  check("sparse_cells_not_reported_as_percentages", all(is.na(cells$reportable_percentage[sparse])))
  check("displayable_percentages", all(abs(cells$reportable_percentage[!sparse] -
         100 * cells$numerator_n[!sparse] / cells$denominator_n[!sparse]) < 1e-8))
}
table <- do.call(rbind, checks)
utils::write.table(table, file.path(qa, "checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
status <- if (all(table$pass)) "DESCRIPTIVE_INDEPENDENT_QA_PASS" else "DESCRIPTIVE_INDEPENDENT_QA_FAIL"
writeLines(status, file.path(qa, "status.txt"))
cat(status, "checks=", nrow(table), "\n")
if (!all(table$pass)) quit(save = "no", status = 1L)
