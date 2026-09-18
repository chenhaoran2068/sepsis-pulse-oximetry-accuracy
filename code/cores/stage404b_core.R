s404b_reverse_pattern <- function(x, threshold, bootstrap_replicates, seed) {
  x <- data.table::as.data.table(x)
  s404_stop_if(!threshold %in% c(88, 92), "STAGE404B_THRESHOLD_INVALID")
  denominator_pairs <- x[spo2_saturation_percent >= threshold]
  numerator_pairs <- denominator_pairs[sao2_saturation_percent < 88]
  if (!nrow(denominator_pairs)) {
    return(data.table::data.table(
      displayed_spo2_threshold_ge = threshold,
      accepted_pair_n = nrow(x),
      conditional_denominator_spo2_ge_n = 0L,
      conditional_numerator_sao2_lt88_n = 0L,
      conditional_proportion = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      denominator_patient_n = 0L, numerator_patient_n = 0L,
      bootstrap_requested = 0L, bootstrap_successful = 0L,
      bootstrap_failed = 0L, bootstrap_seed = NA_integer_,
      cell_status = "NO_DENOMINATOR"
    ))
  }
  by_patient <- denominator_pairs[, .(
    numerator = sum(sao2_saturation_percent < 88),
    denominator = .N
  ), by = patient_key_internal]
  boot <- s404_cluster_proportion(by_patient, bootstrap_replicates, seed)
  data.table::data.table(
    displayed_spo2_threshold_ge = threshold,
    accepted_pair_n = nrow(x),
    conditional_denominator_spo2_ge_n = nrow(denominator_pairs),
    conditional_numerator_sao2_lt88_n = nrow(numerator_pairs),
    conditional_proportion = nrow(numerator_pairs) / nrow(denominator_pairs),
    ci_lower = boot$ci_lower, ci_upper = boot$ci_upper,
    denominator_patient_n = data.table::uniqueN(denominator_pairs$patient_key_internal),
    numerator_patient_n = data.table::uniqueN(numerator_pairs$patient_key_internal),
    bootstrap_requested = boot$requested,
    bootstrap_successful = boot$successful,
    bootstrap_failed = boot$failed,
    bootstrap_seed = boot$seed,
    cell_status = "ESTIMATED"
  )
}

s404b_forbidden_aggregate_names <- function(names_vector) {
  any(grepl(
    "patient_key|stay_key|event_ordinal|timestamp|relative_icu_minutes|saturation_percent",
    names_vector, ignore.case = TRUE
  ))
}
