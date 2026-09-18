mi1_stop <- function(code) stop(code, call. = FALSE)

mi1_require <- function(x, fields, code) {
  if (!is.data.frame(x) || !all(fields %in% names(x))) mi1_stop(code)
}

mi1_pair_key <- function(x) {
  do.call(paste, c(x[ci1_pair_keys()], sep = "\r"))
}

mi1_lmm_frame <- function(analytic_pairs, matched_pairs, cohort, window_minutes) {
  mi1_require(analytic_pairs, c(ci1_pair_keys(), "cohort", "pair_window_minutes",
    "bias_spo2_minus_sao2", "paired_mean_saturation_c90", "age_years_numeric",
    "age_interval", "sex_class", "cci_value", "sofa_total",
    paste0(ci1_metric_names(), "_candidate_value")), "MI1_ANALYTIC_INTERFACE_INVALID")
  mi1_require(matched_pairs, c(ci1_pair_keys(), "cohort", "pair_window_minutes",
    "sao2_relative_icu_minutes"), "MI1_PAIR_TIME_INTERFACE_INVALID")
  if (any(analytic_pairs$cohort != cohort) || any(matched_pairs$cohort != cohort) ||
      any(analytic_pairs$pair_window_minutes != window_minutes) ||
      any(matched_pairs$pair_window_minutes != window_minutes) ||
      anyDuplicated(mi1_pair_key(analytic_pairs)) ||
      anyDuplicated(mi1_pair_key(matched_pairs)) ||
      !setequal(mi1_pair_key(analytic_pairs), mi1_pair_key(matched_pairs)))
    mi1_stop("MI1_PAIR_SCOPE_OR_ALIGNMENT_INVALID")
  index <- match(mi1_pair_key(analytic_pairs), mi1_pair_key(matched_pairs))
  minutes <- matched_pairs$sao2_relative_icu_minutes[index]
  if (anyNA(minutes) || any(!is.finite(minutes)) || any(minutes < 0))
    mi1_stop("MI1_SAO2_ANCHOR_TIME_INVALID")
  out <- data.table::as.data.table(analytic_pairs)
  out[, sao2_anchor_icu_hours := as.numeric(minutes) / 60]
  out <- s405m_derive_age_fields(out, cohort)
  required <- c("patient_key_internal", "stay_key_internal", "bias_spo2_minus_sao2",
    "paired_mean_saturation_c90", "sao2_anchor_icu_hours",
    s405m_stay_fields(cohort), s405m_pair_fields())
  if (!all(required %in% names(out))) mi1_stop("MI1_LMM_FIELDS_MISSING")
  if (anyNA(out$bias_spo2_minus_sao2) || anyNA(out$paired_mean_saturation_c90))
    mi1_stop("MI1_STRUCTURAL_FIELDS_MISSING")
  out
}

mi1_rf_candidate <- function(analytic_pairs, cohort) {
  mi1_require(analytic_pairs, c("cohort", "pair_window_minutes", "patient_key_internal",
    "stay_key_internal", "pair_key_internal", "bias_spo2_minus_sao2",
    "spo2_saturation_percent", "sao2_saturation_percent", "age_years_numeric",
    "age_interval", "age_model_representation", "sex_class", "cci_value",
    "sofa_total", "sofa_status", ci1_metric_columns(),
    setdiff(ci1_pair_columns(), c(ci1_pair_keys(), "pair_key_internal", "sofa_total",
      "sofa_status", ci1_metric_columns()))), "MI1_RF_INTERFACE_INVALID")
  if (any(analytic_pairs$cohort != cohort) ||
      any(analytic_pairs$pair_window_minutes != 60L))
    mi1_stop("MI1_RF_REQUIRES_60_MINUTE_PAIRS")
  out <- data.frame(
    cohort = analytic_pairs$cohort,
    analysis_patient_key = analytic_pairs$patient_key_internal,
    analysis_stay_key = analytic_pairs$stay_key_internal,
    pair_uid = analytic_pairs$pair_key_internal,
    bias_spo2_minus_sao2 = analytic_pairs$bias_spo2_minus_sao2,
    spo2_percent = analytic_pairs$spo2_saturation_percent,
    sao2_percent = analytic_pairs$sao2_saturation_percent,
    age_years_numeric = analytic_pairs$age_years_numeric,
    age_interval = analytic_pairs$age_interval,
    age_model_representation = analytic_pairs$age_model_representation,
    sex_class = analytic_pairs$sex_class,
    cci_value = analytic_pairs$cci_value,
    sofa_total = analytic_pairs$sofa_total,
    sofa_status = analytic_pairs$sofa_status,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  for (name in ci1_metric_names()) {
    for (suffix in c("candidate_value", "canonical_unit", "linkage_class", "value_state")) {
      field <- paste0(name, "_", suffix)
      out[[field]] <- analytic_pairs[[field]]
    }
  }
  for (name in ci1_treatment_names()) {
    for (suffix in c("state", "scope")) {
      field <- paste0(name, "_", suffix)
      out[[field]] <- analytic_pairs[[field]]
    }
  }
  out <- out[r211b_feature_role_manifest()$field]
  r211b_validate_candidate(out, cohort)
  out
}
