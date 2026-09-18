# Private synthetic-only adapter. The R211B core must already be sourced.
private_build_model_inputs <- function(pairs) {
  stopifnot(is.data.frame(pairs), nrow(pairs) > 0L, all(pairs$cohort == "MIMIC"),
    all(pairs$pair_window_minutes == 60))
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  pairs$pair_key_internal <- sprintf("synthetic_model_pair_%05d", seq_len(nrow(pairs)))
  stay <- unique(pairs[c("patient_key_internal", "stay_key_internal")])
  stay$synthetic_index <- seq_len(nrow(stay))
  topcoded <- stay$synthetic_index %% 20L == 0L
  stay$age_years_numeric <- ifelse(topcoded, NA_real_, 35 + (stay$synthetic_index * 7L) %% 54L)
  stay$age_interval <- ifelse(topcoded, "90+", NA_character_)
  stay$age_model_representation <- ifelse(topcoded, "source_interval", "numeric")
  stay$sex_class <- ifelse(stay$synthetic_index %% 2L == 0L, "female", "male")
  stay$cci_value <- as.numeric((stay$synthetic_index * 3L) %% 8L)
  core <- stay[c("patient_key_internal", "stay_key_internal", "age_years_numeric",
    "age_interval", "age_model_representation", "sex_class", "cci_value")]
  index <- match(pairs$stay_key_internal, stay$stay_key_internal)
  j <- as.integer(pairs$sao2_stable_event_ordinal)
  stopifnot(!anyNA(index), all(j %in% 1:4))
  sofa <- data.frame(stay_key_internal = pairs$stay_key_internal,
    hour_end_min = pairs$sao2_relative_icu_minutes,
    sofa_total = as.numeric(pmax(0, round(8 + 4 * sin(index * 1.713) + 2 * cos(index * j * 0.477)))),
    sofa_status = "available", stringsAsFactors = FALSE)
  metric <- pairs[c("patient_key_internal", "stay_key_internal", "sao2_stable_event_ordinal")]
  metric$in_m60 <- TRUE
  values <- list(
    pao2 = as.numeric(60 + (index * 7L + j * 5L) %% 70L),
    paco2 = as.numeric(42 + 8 * sin(index * 1.337) + 2 * cos(index * j * 0.827)),
    ph = as.numeric(7.22 + ((index * 9L + j * 4L) %% 32L) / 100),
    lactate = as.numeric(0.7 + ((index * 5L + j * 2L) %% 42L) / 10),
    hb = as.numeric(8 + ((index * 7L + j * 3L) %% 73L) / 10)
  )
  for (name in names(values)) {
    metric[[paste0(name, "_candidate_value")]] <- values[[name]]
    metric[[paste0(name, "_canonical_unit")]] <- r211b_expected_metric_units()[[name]]
    metric[[paste0(name, "_linkage_class")]] <- "same_specimen"
    metric[[paste0(name, "_value_state")]] <- "eligible_pending_candidate_build"
  }
  clinical <- metric[c("patient_key_internal", "stay_key_internal", "sao2_stable_event_ordinal", "in_m60")]
  clinical$concurrent_vasoactive_medication_use <- ifelse(index %% 3L == 0L, "yes", "unknown")
  clinical$current_invasive_mechanical_ventilation <- ifelse(index %% 4L == 0L, "yes", "unknown")
  base <- r211b_build_candidate("MIMIC", pairs, core, sofa, metric, clinical, TRUE, TRUE)
  base$paired_mean_saturation_c90 <- (base$spo2_percent + base$sao2_percent) / 2 - 90
  base$age_value_per_10y <- ifelse(is.na(base$age_years_numeric), 9, base$age_years_numeric / 10)
  base$age_topcoded_indicator <- as.integer(!is.na(base$age_interval))
  list(pairs = pairs, core = core, candidate = base)
}
