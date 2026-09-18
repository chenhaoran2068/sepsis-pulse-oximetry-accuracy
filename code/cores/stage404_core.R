s404_cohorts <- function() c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")

s404_display_names <- function() c(
  MIMIC = "MIMIC", Amsterdam = "AmsterdamUMCdb", eICU = "eICU",
  SICDB = "SICdb", Lianyungang = "Lianyungang"
)

s404_strata <- function() c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")

s404_stop_if <- function(condition, code) {
  if (isTRUE(condition)) stop(code, call. = FALSE)
  invisible(TRUE)
}

s404_sha256 <- function(path) unname(tools::sha256sum(path))

s404_validate_output_root <- function(path) {
  s404_stop_if(file.exists(path) || dir.exists(path), "STAGE404_OUTPUT_ROOT_ALREADY_EXISTS")
  invisible(TRUE)
}

s404_spo2_stratum <- function(spo2) {
  ifelse(
    spo2 < 88, "70%-<88%",
    ifelse(spo2 < 92, "88%-<92%", ifelse(spo2 <= 96, "92%-<=96%", ">96%-100%")))
}

s404_source_stratum_code <- function(spo2) {
  ifelse(
    spo2 < 88, "70_to_lt88",
    ifelse(spo2 < 92, "88_to_lt92", ifelse(spo2 <= 96, "92_to_96", "gt96_to_100"))
  )
}

s404_source_to_display_stratum <- function(source_code) {
  mapping <- c(
    "70_to_lt88" = "70%-<88%",
    "88_to_lt92" = "88%-<92%",
    "92_to_96" = "92%-<=96%",
    "gt96_to_100" = ">96%-100%"
  )
  unname(mapping[as.character(source_code)])
}

s404_validate_full_pairs <- function(x, cohort, window) {
  required <- c(
    "pair_key_internal", "cohort", "patient_key_internal", "stay_key_internal",
    "spo2_saturation_percent", "sao2_saturation_percent", "absolute_lag_minutes",
    "pair_window_minutes", "contains_possible_transient_artifact"
  )
  s404_stop_if(length(setdiff(required, names(x))) > 0L, "STAGE404_FULL_PAIR_COLUMNS_MISSING")
  s404_stop_if(!nrow(x), "STAGE404_FULL_PAIR_INPUT_EMPTY")
  s404_stop_if(anyNA(x[, ..required]), "STAGE404_FULL_PAIR_REQUIRED_VALUE_MISSING")
  s404_stop_if(data.table::uniqueN(x$pair_key_internal) != nrow(x), "STAGE404_FULL_PAIR_KEY_NOT_UNIQUE")
  s404_stop_if(data.table::uniqueN(x$cohort) != 1L || as.character(x$cohort[[1L]]) != cohort, "STAGE404_FULL_PAIR_COHORT_MISMATCH")
  allowed_lag <- if (window == "M60") 60 else 5
  s404_stop_if(any(as.numeric(x$pair_window_minutes) != allowed_lag), "STAGE404_FULL_PAIR_WINDOW_MISMATCH")
  s404_stop_if(any(as.numeric(x$absolute_lag_minutes) < 0 | as.numeric(x$absolute_lag_minutes) > allowed_lag + 1e-8), "STAGE404_FULL_PAIR_LAG_INVALID")
  s404_stop_if(any(x$spo2_saturation_percent < 70 | x$spo2_saturation_percent > 100), "STAGE404_FULL_SPO2_RANGE_INVALID")
  s404_stop_if(any(x$sao2_saturation_percent < 70 | x$sao2_saturation_percent > 100), "STAGE404_FULL_SAO2_RANGE_INVALID")
  s404_stop_if(any(as.logical(x$contains_possible_transient_artifact)), "STAGE404_FINAL_PAIR_ARTIFACT_PRESENT")
  invisible(TRUE)
}

s404_validate_q21_pairs <- function(x, cohort, window) {
  required <- c(
    "cohort", "patient_key_internal", "stay_key_internal", "spo2_stable_event_ordinal",
    "sao2_stable_event_ordinal", "spo2_relative_icu_minutes", "sao2_relative_icu_minutes",
    "spo2_saturation_percent", "sao2_saturation_percent", "absolute_lag_minutes",
    "pair_window_minutes", "contains_possible_transient_artifact", "relative_icu_day",
    "spo2_display_stratum", "sao2_lt88"
  )
  s404_stop_if(length(setdiff(required, names(x))) > 0L, "STAGE404_Q21_PAIR_COLUMNS_MISSING")
  s404_stop_if(!nrow(x), "STAGE404_Q21_PAIR_INPUT_EMPTY")
  s404_stop_if(anyNA(x[, ..required]), "STAGE404_Q21_PAIR_REQUIRED_VALUE_MISSING")
  s404_stop_if(data.table::uniqueN(x, by = c("patient_key_internal", "stay_key_internal", "spo2_stable_event_ordinal", "sao2_stable_event_ordinal")) != nrow(x), "STAGE404_Q21_PAIR_IDENTITY_NOT_UNIQUE")
  s404_stop_if(data.table::uniqueN(x$cohort) != 1L || as.character(x$cohort[[1L]]) != cohort, "STAGE404_Q21_COHORT_MISMATCH")
  allowed_lag <- if (window == "M60") 60 else 5
  s404_stop_if(any(as.numeric(x$pair_window_minutes) != allowed_lag), "STAGE404_Q21_WINDOW_MISMATCH")
  s404_stop_if(any(x$spo2_saturation_percent < 70 | x$spo2_saturation_percent > 100), "STAGE404_Q21_SPO2_RANGE_INVALID")
  s404_stop_if(any(x$sao2_saturation_percent < 70 | x$sao2_saturation_percent > 100), "STAGE404_Q21_SAO2_RANGE_INVALID")
  s404_stop_if(any(as.logical(x$contains_possible_transient_artifact)), "STAGE404_Q21_FINAL_PAIR_ARTIFACT_PRESENT")
  s404_stop_if(any(as.numeric(x$spo2_relative_icu_minutes) < 0 | as.numeric(x$sao2_relative_icu_minutes) < 0 | as.numeric(x$spo2_relative_icu_minutes) > 10080 | as.numeric(x$sao2_relative_icu_minutes) > 10080), "STAGE404_Q21_SCOPE_OUTSIDE_0_168H")
  s404_stop_if(any(as.integer(x$relative_icu_day) < 1L | as.integer(x$relative_icu_day) > 7L), "STAGE404_Q21_DAY_INVALID")
  calculated_source_stratum <- s404_source_stratum_code(as.numeric(x$spo2_saturation_percent))
  s404_stop_if(!all(calculated_source_stratum == as.character(x$spo2_display_stratum)), "STAGE404_Q21_STRATUM_MISMATCH")
  s404_stop_if(!all(as.logical(x$sao2_lt88) == (as.numeric(x$sao2_saturation_percent) < 88)), "STAGE404_Q21_SAO2_FLAG_MISMATCH")
  invisible(TRUE)
}

s404_cluster_proportion <- function(patient_summary, bootstrap_replicates, seed) {
  patient_summary <- data.table::as.data.table(patient_summary)
  s404_stop_if(!all(c("numerator", "denominator") %in% names(patient_summary)), "STAGE404_BOOTSTRAP_INPUT_COLUMNS_MISSING")
  s404_stop_if(nrow(patient_summary) < 1L || sum(patient_summary$denominator) < 1L, "STAGE404_BOOTSTRAP_EMPTY_DENOMINATOR")
  s404_stop_if(any(patient_summary$numerator < 0 | patient_summary$denominator < 0 | patient_summary$numerator > patient_summary$denominator), "STAGE404_BOOTSTRAP_SUBSET_INVALID")
  numerator <- as.numeric(patient_summary$numerator)
  denominator <- as.numeric(patient_summary$denominator)
  patient_n <- length(numerator)
  estimate <- sum(numerator) / sum(denominator)
  values <- numeric(bootstrap_replicates)
  valid <- logical(bootstrap_replicates)
  set.seed(as.integer(seed))
  cursor <- 1L
  while (cursor <= bootstrap_replicates) {
    block_n <- min(250L, bootstrap_replicates - cursor + 1L)
    frequencies <- stats::rmultinom(block_n, size = patient_n, prob = rep(1 / patient_n, patient_n))
    boot_num <- as.numeric(crossprod(numerator, frequencies))
    boot_den <- as.numeric(crossprod(denominator, frequencies))
    index <- cursor:(cursor + block_n - 1L)
    valid[index] <- boot_den > 0
    values[index[valid[index]]] <- boot_num[valid[index]] / boot_den[valid[index]]
    cursor <- cursor + block_n
  }
  successful <- sum(valid)
  list(
    estimate = estimate,
    ci_lower = if (successful) as.numeric(stats::quantile(values[valid], 0.025, names = FALSE, type = 7)) else NA_real_,
    ci_upper = if (successful) as.numeric(stats::quantile(values[valid], 0.975, names = FALSE, type = 7)) else NA_real_,
    requested = bootstrap_replicates,
    successful = successful,
    failed = bootstrap_replicates - successful,
    seed = as.integer(seed)
  )
}

s404_overall_pattern <- function(x, threshold, bootstrap_replicates, seed) {
  x <- data.table::as.data.table(x)
  accepted_patient_n <- data.table::uniqueN(x$patient_key_internal)
  accepted_stay_n <- data.table::uniqueN(x$stay_key_internal)
  at_risk <- x[sao2_saturation_percent < 88]
  event <- at_risk[spo2_saturation_percent >= threshold]
  by_patient <- at_risk[, .(
    numerator = sum(spo2_saturation_percent >= threshold),
    denominator = .N
  ), by = patient_key_internal]
  boot <- if (nrow(at_risk)) s404_cluster_proportion(by_patient, bootstrap_replicates, seed) else NULL
  data.table::data.table(
    threshold_spo2_ge = threshold,
    accepted_pair_n = nrow(x),
    pair_denominator_sao2_lt88_n = nrow(at_risk),
    pair_numerator_n = nrow(event),
    pair_conditional_proportion = if (nrow(at_risk)) nrow(event) / nrow(at_risk) else NA_real_,
    pair_ci_lower = if (is.null(boot)) NA_real_ else boot$ci_lower,
    pair_ci_upper = if (is.null(boot)) NA_real_ else boot$ci_upper,
    accepted_patient_n = accepted_patient_n,
    patient_numerator_at_least_one_episode_n = data.table::uniqueN(event$patient_key_internal),
    accepted_stay_n = accepted_stay_n,
    stay_numerator_at_least_one_episode_n = data.table::uniqueN(event$stay_key_internal),
    bootstrap_requested = if (is.null(boot)) 0L else boot$requested,
    bootstrap_successful = if (is.null(boot)) 0L else boot$successful,
    bootstrap_failed = if (is.null(boot)) 0L else boot$failed,
    bootstrap_seed = if (is.null(boot)) NA_integer_ else boot$seed,
    cell_status = if (!nrow(at_risk)) "NO_DENOMINATOR" else "ESTIMATED"
  )
}

s404_day_stratum_cell <- function(x, bootstrap_replicates, seed) {
  x <- data.table::as.data.table(x)
  numerator <- sum(as.logical(x$sao2_lt88))
  denominator <- nrow(x)
  patient_n <- data.table::uniqueN(x$patient_key_internal)
  stay_n <- data.table::uniqueN(x$stay_key_internal)
  if (!denominator) {
    return(data.table::data.table(
      numerator_n = 0L, denominator_n = 0L, patient_n = 0L, stay_n = 0L,
      conditional_proportion = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      bootstrap_requested = 0L, bootstrap_successful = 0L, bootstrap_failed = 0L,
      bootstrap_seed = NA_integer_, cell_status = "NO_DENOMINATOR"
    ))
  }
  if (denominator < 20L || patient_n < 10L) {
    return(data.table::data.table(
      numerator_n = numerator, denominator_n = denominator, patient_n = patient_n, stay_n = stay_n,
      conditional_proportion = numerator / denominator, ci_lower = NA_real_, ci_upper = NA_real_,
      bootstrap_requested = 0L, bootstrap_successful = 0L, bootstrap_failed = 0L,
      bootstrap_seed = NA_integer_, cell_status = "SPARSE_NOT_PLOTTED"
    ))
  }
  by_patient <- x[, .(numerator = sum(as.logical(sao2_lt88)), denominator = .N), by = patient_key_internal]
  boot <- s404_cluster_proportion(by_patient, bootstrap_replicates, seed)
  data.table::data.table(
    numerator_n = numerator, denominator_n = denominator, patient_n = patient_n, stay_n = stay_n,
    conditional_proportion = numerator / denominator, ci_lower = boot$ci_lower, ci_upper = boot$ci_upper,
    bootstrap_requested = boot$requested, bootstrap_successful = boot$successful, bootstrap_failed = boot$failed,
    bootstrap_seed = boot$seed, cell_status = "DISPLAYABLE"
  )
}

s404_forbidden_aggregate_names <- function(names_vector) {
  any(grepl("patient_key|stay_key|event_ordinal|relative_icu_minutes|timestamp|saturation_percent", names_vector, ignore.case = TRUE))
}
