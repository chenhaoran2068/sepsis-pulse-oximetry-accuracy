s404c_pattern_summary <- function(x, threshold) {
  x <- data.table::as.data.table(x)
  s404_stop_if(!threshold %in% c(88, 92), "STAGE404C_THRESHOLD_INVALID")
  x[, event := sao2_saturation_percent < 88 & spo2_saturation_percent >= threshold]
  patient_counts <- x[, .(event_pair_n = sum(event)), by = patient_key_internal]
  stay_counts <- x[, .(event_pair_n = sum(event)), by = stay_key_internal]
  affected <- patient_counts[event_pair_n > 0]
  affected_stays <- stay_counts[event_pair_n > 0]
  event_total <- sum(patient_counts$event_pair_n)
  patient_summary <- data.table::data.table(
    displayed_spo2_threshold_ge = threshold,
    accepted_patient_n = nrow(patient_counts),
    affected_patient_n = nrow(affected),
    affected_patient_proportion = nrow(affected) / nrow(patient_counts),
    accepted_stay_n = nrow(stay_counts),
    affected_stay_n = nrow(affected_stays),
    affected_stay_proportion = nrow(affected_stays) / nrow(stay_counts),
    observed_event_pair_n = event_total
  )
  if (!nrow(affected)) {
    concentration <- data.table::data.table(
      displayed_spo2_threshold_ge = threshold,
      affected_patient_n = 0L, observed_event_pair_n = 0L,
      event_pairs_per_affected_patient_q1 = NA_real_,
      event_pairs_per_affected_patient_median = NA_real_,
      event_pairs_per_affected_patient_q3 = NA_real_,
      event_pairs_per_affected_patient_max = NA_integer_,
      affected_patients_1_pair_n = 0L, affected_patients_1_pair_proportion = NA_real_,
      affected_patients_2_pairs_n = 0L, affected_patients_2_pairs_proportion = NA_real_,
      affected_patients_ge3_pairs_n = 0L, affected_patients_ge3_pairs_proportion = NA_real_,
      top_10_percent_affected_patient_n = 0L, top_10_percent_event_pair_share = NA_real_,
      cell_status = "NO_AFFECTED_PATIENTS"
    )
    return(list(patient_summary = patient_summary, concentration = concentration))
  }
  counts <- sort(as.integer(affected$event_pair_n), decreasing = TRUE)
  q <- as.numeric(stats::quantile(counts, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 7))
  top_n <- max(1L, ceiling(0.10 * length(counts)))
  concentration <- data.table::data.table(
    displayed_spo2_threshold_ge = threshold,
    affected_patient_n = length(counts), observed_event_pair_n = sum(counts),
    event_pairs_per_affected_patient_q1 = q[[1L]],
    event_pairs_per_affected_patient_median = q[[2L]],
    event_pairs_per_affected_patient_q3 = q[[3L]],
    event_pairs_per_affected_patient_max = max(counts),
    affected_patients_1_pair_n = sum(counts == 1L),
    affected_patients_1_pair_proportion = mean(counts == 1L),
    affected_patients_2_pairs_n = sum(counts == 2L),
    affected_patients_2_pairs_proportion = mean(counts == 2L),
    affected_patients_ge3_pairs_n = sum(counts >= 3L),
    affected_patients_ge3_pairs_proportion = mean(counts >= 3L),
    top_10_percent_affected_patient_n = top_n,
    top_10_percent_event_pair_share = sum(counts[seq_len(top_n)]) / sum(counts),
    cell_status = "AUDITED"
  )
  list(patient_summary = patient_summary, concentration = concentration)
}

s404c_forbidden_aggregate_names <- function(names_vector) {
  any(grepl(
    "patient_key|stay_key|event_ordinal|timestamp|relative_icu_minutes|saturation_percent",
    names_vector, ignore.case = TRUE
  ))
}
