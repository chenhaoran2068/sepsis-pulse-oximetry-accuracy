s405m_cohorts <- function() c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
s405m_windows <- function() c("M60", "M5")
s405m_metric_names <- function() c("pao2", "paco2", "ph", "lactate", "hb")
s405m_key_cols <- function() c("cohort", "patient_key_internal", "stay_key_internal", "sao2_stable_event_ordinal")

s405m_stop_if <- function(condition, code) {
  if (isTRUE(condition)) stop(code, call. = FALSE)
  invisible(TRUE)
}

mice.impute.2l.pmm.scaled <- function(y, ry, x, type, ...) {
  x <- as.matrix(x)
  scalable <- which(type != -2L & type != 0L)
  for (column in scalable) {
    value <- suppressWarnings(as.numeric(x[, column]))
    finite <- is.finite(value)
    observed_unique_n <- length(unique(value[finite]))
    scale_value <- if (sum(finite) > 1L) stats::sd(value[finite]) else NA_real_
    if (observed_unique_n > 2L && is.finite(scale_value) && scale_value > 0) {
      x[, column] <- (value - mean(value[finite])) / scale_value
    }
  }
  miceadds::mice.impute.2l.pmm(y = y, ry = ry, x = x, type = type, ...)
}

s405m_require <- function(x, fields, code) {
  s405m_stop_if(!all(fields %in% names(x)), code)
}

s405m_assert_unique <- function(x, by, code) {
  s405m_stop_if(anyDuplicated(x[, ..by]) > 0L, code)
}

s405m_stay_fields <- function(cohort) {
  if (cohort %in% c("Amsterdam", "SICDB")) {
    c("age_interval", "sex_class", "cci_value")
  } else if (cohort %in% c("MIMIC", "eICU")) {
    c("age_value_per_10y", "age_topcoded_indicator", "sex_class", "cci_value")
  } else {
    c("age_value_per_10y", "sex_class", "cci_value")
  }
}

s405m_age_interval_levels <- function(cohort) {
  switch(cohort,
    Amsterdam = c("18-39", "40-49", "50-59", "60-69", "70-79", "80+"),
    SICDB = c(
      paste0(seq(20, 85, by = 5), "_rounded_5y_bin"),
      "90+_rounded_or_topcoded"
    ),
    stop("S405M_AGE_INTERVAL_LEVELS_NOT_APPLICABLE", call. = FALSE)
  )
}

s405m_normalize_sex <- function(value) {
  sex <- tolower(trimws(as.character(value)))
  sex[is.na(value) | sex %in% c("", "unknown", "unknown_or_missing", "missing", "na", "n/a")] <- NA_character_
  s405m_stop_if(any(!is.na(sex) & !sex %in% c("female", "male")), "S405M_SEX_LEVEL_INVALID")
  sex
}

s405m_prepare_analytic_states <- function(x, cohort) {
  out <- data.table::copy(data.table::as.data.table(x))
  s405m_require(out, c("patient_key_internal", "stay_key_internal", "sex_class"), "S405M_SEX_RECONCILIATION_INTERFACE_MISSING")
  out[, .source_row_order := .I]
  out[, .sex_normalized := s405m_normalize_sex(sex_class)]
  stay <- out[, {
    observed <- unique(.sex_normalized[!is.na(.sex_normalized)])
    s405m_stop_if(length(observed) > 1L, "S405M_SEX_WITHIN_STAY_CONFLICT")
    .(
      source_sex = if (length(observed)) observed[[1L]] else NA_character_,
      source_observed = length(observed) == 1L
    )
  }, by = .(patient_key_internal, stay_key_internal)]
  patient <- stay[, {
    observed <- unique(source_sex[!is.na(source_sex)])
    s405m_stop_if(length(observed) > 1L, "S405M_SEX_WITHIN_PATIENT_CONFLICT")
    .(
      patient_observed_sex = if (length(observed)) observed[[1L]] else NA_character_,
      patient_stay_n = .N
    )
  }, by = patient_key_internal]
  stay <- merge(stay, patient, by = "patient_key_internal", all.x = TRUE, sort = FALSE)
  stay[, reconciled_from_patient := is.na(source_sex) & !is.na(patient_observed_sex)]
  stay[, analytic_sex := data.table::fifelse(reconciled_from_patient, patient_observed_sex, source_sex)]
  s405m_stop_if(
    any(is.na(stay$analytic_sex) & stay$patient_stay_n > 1L),
    "S405M_SEX_ALL_MISSING_MULTI_STAY_REQUIRES_PATIENT_LEVEL_METHOD"
  )
  stay[, analytic_observed := !is.na(analytic_sex)]
  out[stay, on = .(patient_key_internal, stay_key_internal), `:=`(
    sex_class = i.analytic_sex,
    .sex_source_observed = i.source_observed,
    .sex_reconciled_from_patient = i.reconciled_from_patient
  )]
  data.table::setorder(out, .source_row_order)
  out[, c(".source_row_order", ".sex_normalized") := NULL]
  list(
    data = out,
    sex_ledger = stay[, .(
      patient_key_internal, stay_key_internal,
      source_observed, reconciled_from_patient, analytic_observed,
      residual_missing = is.na(analytic_sex)
    )]
  )
}

s405m_derive_age_fields <- function(stays, cohort) {
  s405m_require(stays, c("age_years_numeric", "age_interval"), "S405M_RAW_AGE_INTERFACE_MISSING")
  numeric_age <- suppressWarnings(as.numeric(stays$age_years_numeric))
  interval <- trimws(as.character(stays$age_interval))
  interval[is.na(stays$age_interval) | interval == ""] <- NA_character_

  if (cohort %in% c("MIMIC", "eICU")) {
    topcoded <- !is.na(interval) & grepl("^[0-9]+\\+$", interval)
    s405m_stop_if(any(!is.na(interval) & !topcoded), "S405M_TOPCODE_LABEL_INVALID")
    lower_bound <- rep(NA_real_, length(interval))
    lower_bound[topcoded] <- as.numeric(sub("\\+$", "", interval[topcoded]))
    s405m_stop_if(any(!is.na(numeric_age) & topcoded), "S405M_NUMERIC_AND_TOPCODE_BOTH_PRESENT")
    if (cohort == "MIMIC") {
      s405m_stop_if(any(topcoded & !lower_bound %in% c(90, 91, 92)), "S405M_MIMIC_TOPCODE_LOWER_BOUND_INVALID")
    } else {
      s405m_stop_if(any(topcoded & lower_bound != 90), "S405M_EICU_TOPCODE_LOWER_BOUND_INVALID")
    }
    age_value <- ifelse(!is.na(numeric_age), numeric_age, lower_bound)
    stays[, age_value_per_10y := age_value / 10]
    stays[, age_topcoded_indicator := as.integer(topcoded)]
    s405m_stop_if(anyNA(stays$age_topcoded_indicator), "S405M_TOPCODE_INDICATOR_MISSING")
  } else if (cohort == "Lianyungang") {
    stays[, age_value_per_10y := numeric_age / 10]
  } else if (cohort %in% c("Amsterdam", "SICDB")) {
    s405m_stop_if(any(!is.na(numeric_age)), "S405M_INTERVAL_COHORT_NUMERIC_AGE_PRESENT")
    stays[, age_interval := interval]
  } else {
    stop("S405M_UNKNOWN_COHORT", call. = FALSE)
  }
  stays
}

s405m_pair_fields <- function() {
  c("sofa_total", paste0(s405m_metric_names(), "_candidate_value"))
}

s405m_allowed_predictor_fields <- function(cohort) {
  common <- c(
    s405m_stay_fields(cohort), "sofa_total",
    "pao2_candidate_value", "paco2_candidate_value", "ph_candidate_value"
  )
  additional <- switch(cohort,
    Amsterdam = "hb_candidate_value",
    SICDB = c("lactate_candidate_value", "hb_candidate_value"),
    Lianyungang = c("lactate_candidate_value", "hb_candidate_value"),
    MIMIC = character(), eICU = character(),
    stop("S405M_UNKNOWN_COHORT", call. = FALSE)
  )
  unique(c(common, additional))
}

s405m_metric_sidecar_identity <- function(project_root, cohort) {
  base <- file.path(project_root, "06_data", "03_analysis_ready_datasets", "controlled_candidates")
  repaired_root <- file.path(base, "R212D-PH-VALUE-RANGE-REPAIR-RENAL-SOFA-P1-CANDIDATE-REAL-20260905-229")
  expected_sha256 <- switch(cohort,
    MIMIC = "140C076CE982E6F51041A74ED521EB54D4DF34BB74B66B0950A13B3C252FC1EF",
    Amsterdam = "7F954F03445046F3C28CDC32C7A19259A305D20C5644302C37AB85789FD9D6E7",
    eICU = "B230E4955FDDB46E97998B53E3ECB6535F01252417B7CD3A24B9B76D5D08F63F",
    SICDB = "642F86BEFB0DDF41D02C8248D2C29314159DE0BCF8C78119DD8CF7D1E88F7892",
    Lianyungang = "74D3BD83D824EE83022BFC197736C90C230E705E5CA70CA317D42B929EC2D8E2",
    stop("S405M_UNKNOWN_COHORT", call. = FALSE)
  )
  list(
    cohort = cohort,
    path = file.path(repaired_root, cohort, "time_linked_ph_validity_repaired_candidate_NOT_ADOPTED.parquet"),
    expected_sha256 = expected_sha256,
    authority = "renal_sofa_p1_reentry_ph_repair_independent_qa_pass_pending_owner_gate"
  )
}

s405m_metric_sidecar_path <- function(project_root, cohort) {
  s405m_metric_sidecar_identity(project_root, cohort)$path
}

s405m_assert_file_identity <- function(path, expected_sha256,
                                       missing_code = "S405M_APPROVED_INPUT_MISSING",
                                       mismatch_code = "S405M_APPROVED_INPUT_HASH_MISMATCH") {
  s405m_stop_if(!file.exists(path), missing_code)
  actual_sha256 <- digest::digest(path, file = TRUE, algo = "sha256")
  s405m_stop_if(!identical(tolower(actual_sha256), tolower(expected_sha256)), mismatch_code)
  actual_sha256
}

s405m_assert_metric_sidecar_approved <- function(project_root, cohort) {
  identity <- s405m_metric_sidecar_identity(project_root, cohort)
  identity$actual_sha256 <- s405m_assert_file_identity(
    identity$path,
    identity$expected_sha256,
    missing_code = "S405M_APPROVED_METRIC_SIDECAR_MISSING",
    mismatch_code = "S405M_APPROVED_METRIC_SIDECAR_HASH_MISMATCH"
  )
  identity
}

s405m_assert_retained_ph_range <- function(side) {
  s405m_require(side, "ph_candidate_value", "S405M_PH_INTERFACE_MISSING")
  value <- side$ph_candidate_value
  invalid <- !is.na(value) & (!is.finite(value) | value < 6.3 | value > 8.0)
  s405m_stop_if(any(invalid), "S405M_RETAINED_PH_OUTSIDE_APPROVED_RANGE")
  invisible(TRUE)
}

s405m_metric_validity_map <- function() {
  data.table::data.table(
    metric = s405m_metric_names(),
    value_field = paste0(s405m_metric_names(), "_candidate_value"),
    validity_variable = c("pao2_mmhg", "paco2_mmhg", "ph", "lactate_mmol_l", "hb_g_dl")
  )
}

s405m_validity_core <- function(project_root) {
  path <- file.path(
    project_root, "07_analysis", "controlled_pipeline_candidates",
    "STAGE4-VARIABLE-VALIDITY-AMENDMENT-R1-NODATA-20260906-280",
    "R", "stage4_variable_validity_core.R"
  )
  s405m_stop_if(!file.exists(path), "S405M_VALIDITY_CORE_MISSING")
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  s405m_stop_if(!exists("s4v_validate", envir = env, inherits = FALSE), "S405M_VALIDITY_INTERFACE_MISSING")
  list(path = path, validate = get("s4v_validate", envir = env, inherits = FALSE))
}

s405m_apply_metric_validity <- function(side, project_root) {
  out <- data.table::copy(data.table::as.data.table(side))
  mapping <- s405m_metric_validity_map()
  s405m_require(out, mapping$value_field, "S405M_METRIC_VALIDITY_INPUT_MISSING")
  validity <- s405m_validity_core(project_root)
  for (i in seq_len(nrow(mapping))) {
    value_field <- mapping$value_field[[i]]
    checked <- validity$validate(out[[value_field]], mapping$validity_variable[[i]])
    out[[value_field]] <- checked$value
    out[[paste0(value_field, "_validity_state")]] <- checked$value_state
    out[[paste0(value_field, "_diagnostic_review")]] <- checked$diagnostic_review
  }
  attr(out, "validity_core_path") <- validity$path
  out
}

s405m_metric_validity_ledger <- function(x, cohort, window) {
  mapping <- s405m_metric_validity_map()
  data.table::rbindlist(lapply(seq_len(nrow(mapping)), function(i) {
    value_field <- mapping$value_field[[i]]
    state_field <- paste0(value_field, "_validity_state")
    review_field <- paste0(value_field, "_diagnostic_review")
    s405m_require(x, c(value_field, state_field, review_field), "S405M_METRIC_VALIDITY_LEDGER_INTERFACE")
    state <- as.character(x[[state_field]])
    data.table::data.table(
      cohort = cohort,
      analysis_window = window,
      metric = mapping$metric[[i]],
      pair_n = nrow(x),
      available_n = sum(is.finite(x[[value_field]])),
      invalid_n = sum(grepl("^invalid_", state)),
      missing_or_unavailable_n = sum(state == "missing_or_unavailable"),
      diagnostic_review_n = sum(x[[review_field]] != "none")
    )
  }))
}

s405m_clinical_sidecar_path <- function(project_root, cohort) {
  base <- file.path(project_root, "06_data", "03_analysis_ready_datasets", "controlled_candidates")
  switch(cohort,
    MIMIC = file.path(base, "R08L-MIMIC-CLINICAL-STATE-CANDIDATE-R1-REAL-20260827-198", "sidecars", "MIMIC", "concurrent_clinical_state_candidate_NOT_ADOPTED.parquet"),
    Amsterdam = file.path(base, "R08M-AMSTERDAM-CLINICAL-STATE-CANDIDATE-REAL-20260827-200", "sidecars", "Amsterdam", "concurrent_clinical_state_candidate_NOT_ADOPTED.parquet"),
    eICU = file.path(base, "R08N-EICU-CLINICAL-STATE-CANDIDATE-REAL-20260827-202", "sidecars", "eICU", "concurrent_clinical_state_candidate_NOT_ADOPTED.parquet"),
    SICDB = file.path(base, "R08O-SICDB-VASOACTIVE-CANDIDATE-REAL-20260827-204", "sidecars", "SICDB", "concurrent_vasoactive_medication_candidate_NOT_ADOPTED.parquet"),
    Lianyungang = NA_character_,
    stop("S405M_UNKNOWN_COHORT", call. = FALSE)
  )
}

s405m_attach_sofa <- function(pairs, sofas) {
  p <- data.table::copy(pairs)[, anchor_min := sao2_relative_icu_minutes]
  h <- data.table::copy(sofas)[!is.na(hour_end_min) & !is.na(sofa_total)]
  data.table::setkey(h, cohort, patient_key_internal, stay_key_internal, hour_end_min)
  linked <- h[p,
    on = .(cohort, patient_key_internal, stay_key_internal, hour_end_min <= anchor_min),
    mult = "last", nomatch = NA,
    .(pair_uid = i.pair_uid, sofa_total_linked = x.sofa_total,
      sofa_hour_end_min = x.hour_end_min, anchor_min = i.anchor_min)
  ]
  s405m_assert_unique(linked, "pair_uid", "S405M_SOFA_LINK_DUPLICATE")
  linked[, sofa_delta_min := anchor_min - sofa_hour_end_min]
  linked[is.na(sofa_delta_min) | sofa_delta_min < 0 | sofa_delta_min > 120,
    `:=`(sofa_total_linked = NA_real_, sofa_delta_min = NA_real_)
  ]
  merge(
    p[, setdiff(names(p), "anchor_min"), with = FALSE],
    linked[, .(pair_uid, sofa_total = sofa_total_linked, sofa_delta_min)],
    by = "pair_uid", all.x = TRUE, sort = FALSE
  )
}

s405m_attach_metrics <- function(pairs, project_root, cohort, window) {
  identity <- s405m_assert_metric_sidecar_approved(project_root, cohort)
  path <- identity$path
  side <- data.table::as.data.table(arrow::read_parquet(path, as_data_frame = TRUE))
  keys <- s405m_key_cols()
  values <- paste0(s405m_metric_names(), "_candidate_value")
  for (field in c(values, "in_m60", "in_m5")) {
    if (!field %in% names(side)) {
      match_name <- names(side)[tolower(names(side)) == tolower(field)]
      if (length(match_name) == 1L) data.table::setnames(side, match_name, field)
    }
  }
  s405m_require(side, c(keys, values, "in_m60", "in_m5"), "S405M_METRIC_INTERFACE_MISSING")
  side <- s405m_apply_metric_validity(side, project_root)
  s405m_assert_retained_ph_range(side)
  states <- paste0(values, "_validity_state")
  reviews <- paste0(values, "_diagnostic_review")
  flag <- if (window == "M60") "in_m60" else "in_m5"
  side <- side[get(flag) %in% TRUE, c(keys, values, states, reviews), with = FALSE]
  s405m_assert_unique(side, keys, "S405M_METRIC_DUPLICATE")
  out <- merge(pairs, side, by = keys, all.x = TRUE, sort = FALSE)
  for (i in seq_along(values)) {
    out[is.na(get(states[[i]])), (states[[i]]) := "missing_or_unavailable"]
    out[is.na(get(reviews[[i]])), (reviews[[i]]) := "none"]
  }
  s405m_stop_if(nrow(out) != nrow(pairs), "S405M_METRIC_JOIN_CARDINALITY")
  out
}

s405m_attach_clinical_state <- function(pairs, project_root, cohort, window) {
  path <- s405m_clinical_sidecar_path(project_root, cohort)
  if (is.na(path)) {
    pairs[, `:=`(
      concurrent_vasoactive_medication_use_state = NA_character_,
      current_invasive_mechanical_ventilation_state = NA_character_,
      clinical_state_scope = "not_available_for_cohort"
    )]
    return(pairs)
  }
  s405m_stop_if(!file.exists(path), "S405M_CLINICAL_SIDECAR_MISSING")
  side <- data.table::as.data.table(arrow::read_parquet(path, as_data_frame = TRUE))
  keys <- s405m_key_cols()
  flag <- if (window == "M60") "in_m60" else "in_m5"
  clinical <- intersect(c("concurrent_vasoactive_medication_use", "current_invasive_mechanical_ventilation"), names(side))
  s405m_require(side, c(keys, flag), "S405M_CLINICAL_INTERFACE_MISSING")
  side <- side[get(flag) %in% TRUE, unique(c(keys, clinical)), with = FALSE]
  if (length(clinical)) data.table::setnames(side, clinical, paste0(clinical, "_state"))
  s405m_assert_unique(side, keys, "S405M_CLINICAL_DUPLICATE")
  out <- merge(pairs, side, by = keys, all.x = TRUE, sort = FALSE)
  if (!"concurrent_vasoactive_medication_use_state" %in% names(out)) out[, concurrent_vasoactive_medication_use_state := NA_character_]
  if (!"current_invasive_mechanical_ventilation_state" %in% names(out)) out[, current_invasive_mechanical_ventilation_state := NA_character_]
  out[, clinical_state_scope := "yes_unknown_not_core_predictor"]
  s405m_stop_if(nrow(out) != nrow(pairs), "S405M_CLINICAL_JOIN_CARDINALITY")
  out
}

s405m_read_pair_frame <- function(project_root, cohort, window) {
  s405m_stop_if(!cohort %in% s405m_cohorts(), "S405M_UNKNOWN_COHORT")
  s405m_stop_if(!window %in% s405m_windows(), "S405M_UNKNOWN_WINDOW")
  base <- file.path(project_root, "11_qa", "R06C-CORE-FULL-VARIABLE-RENAL-SOFA-P1-R1-CANDIDATE-REAL-20260905-197", "candidate_output", cohort)
  pair_path <- file.path(base, paste0("pairs_", window, "_NOT_ADOPTED.parquet"))
  stay_path <- file.path(base, "stay_core_NOT_ADOPTED.parquet")
  sofa_path <- file.path(base, "hourly_sofa_for_candidate_stays_NOT_ADOPTED.parquet")
  s405m_stop_if(!all(file.exists(c(pair_path, stay_path, sofa_path))), "S405M_BASE_INPUT_MISSING")
  pairs <- data.table::as.data.table(arrow::read_parquet(pair_path, as_data_frame = TRUE))
  stays <- data.table::as.data.table(arrow::read_parquet(stay_path, as_data_frame = TRUE))
  sofas <- data.table::as.data.table(arrow::read_parquet(sofa_path, as_data_frame = TRUE))
  s405m_require(pairs, c("cohort", "patient_key_internal", "stay_key_internal", "spo2_stable_event_ordinal", "sao2_stable_event_ordinal", "spo2_saturation_percent", "sao2_saturation_percent", "sao2_relative_icu_minutes"), "S405M_PAIR_INTERFACE_MISSING")
  s405m_require(stays, c("cohort", "patient_key_internal", "stay_key_internal", "age_years_numeric", "age_interval", "sex_class", "cci_value"), "S405M_STAY_INTERFACE_MISSING")
  s405m_require(sofas, c("cohort", "patient_key_internal", "stay_key_internal", "hour_end_min", "sofa_total"), "S405M_SOFA_INTERFACE_MISSING")
  s405m_assert_unique(stays, c("cohort", "patient_key_internal", "stay_key_internal"), "S405M_STAY_DUPLICATE")
  stays <- s405m_derive_age_fields(stays, cohort)
  pairs[, pair_uid := paste(
    cohort, patient_key_internal, stay_key_internal,
    spo2_stable_event_ordinal, sao2_stable_event_ordinal,
    sep = "::"
  )]
  s405m_assert_unique(pairs, "pair_uid", "S405M_PAIR_UID_DUPLICATE")
  pairs <- merge(
    pairs,
    stays[, c("cohort", "patient_key_internal", "stay_key_internal", s405m_stay_fields(cohort)), with = FALSE],
    by = c("cohort", "patient_key_internal", "stay_key_internal"), all.x = TRUE, sort = FALSE
  )
  s405m_stop_if(nrow(pairs) == 0L, "S405M_EMPTY_PAIR_FRAME")
  pairs[, paired_mean_saturation_c90 := ((spo2_saturation_percent + sao2_saturation_percent) / 2) - 90]
  pairs[, bias_spo2_minus_sao2 := spo2_saturation_percent - sao2_saturation_percent]
  pairs[, sao2_anchor_icu_hours := sao2_relative_icu_minutes / 60]
  s405m_stop_if(anyNA(pairs$bias_spo2_minus_sao2) || anyNA(pairs$paired_mean_saturation_c90) || anyNA(pairs$sao2_anchor_icu_hours), "S405M_OBSERVED_AUXILIARY_MISSING")
  pairs <- s405m_attach_sofa(pairs, sofas)
  pairs <- s405m_attach_metrics(pairs, project_root, cohort, window)
  pairs <- s405m_attach_clinical_state(pairs, project_root, cohort, window)
  pairs
}

s405m_stay_table <- function(x, fields) {
  x[, lapply(.SD, function(v) {
    observed <- unique(v[!is.na(v)])
    s405m_stop_if(length(observed) > 1L, "S405M_STAY_LEVEL_VALUE_CONFLICT")
    if (!length(observed)) v[NA_integer_] else observed[[1L]]
  }), by = .(patient_key_internal, stay_key_internal), .SDcols = fields]
}

s405m_coverage_ledger <- function(x, cohort) {
  prepared <- s405m_prepare_analytic_states(x, cohort)
  x <- prepared$data
  stay_fields <- s405m_stay_fields(cohort)
  pair_fields <- s405m_pair_fields()
  stays <- s405m_stay_table(x, stay_fields)
  stay_rows <- data.table::rbindlist(lapply(stay_fields, function(field) {
    observed <- sum(!is.na(stays[[field]]))
    source_observed <- if (field == "sex_class") sum(prepared$sex_ledger$source_observed) else observed
    reconciled <- if (field == "sex_class") sum(prepared$sex_ledger$reconciled_from_patient) else 0L
    data.table::data.table(
      feature = field, variable_level = "study_icu_stay",
      source_observed_n = source_observed,
      patient_reconciled_n = reconciled,
      observed_n = observed, denominator_n = nrow(stays),
      observed_coverage = observed / nrow(stays)
    )
  }))
  pair_rows <- data.table::rbindlist(lapply(pair_fields, function(field) {
    observed <- sum(!is.na(x[[field]]))
    data.table::data.table(
      feature = field, variable_level = "pair",
      source_observed_n = observed,
      patient_reconciled_n = 0L,
      observed_n = observed, denominator_n = nrow(x),
      observed_coverage = observed / nrow(x)
    )
  }))
  out <- data.table::rbindlist(list(stay_rows, pair_rows))
  out[, coverage_eligible := observed_coverage >= 0.70]
  out[, recipe_allowed := feature %in% s405m_allowed_predictor_fields(cohort)]
  out[, core_model_eligible := coverage_eligible & recipe_allowed]
  out[, omission_reason := data.table::fcase(
    !recipe_allowed, "not_in_prespecified_core_recipe",
    !coverage_eligible, "coverage_below_70_percent_at_variable_level",
    default = ""
  )]
  out
}

s405m_select_predictors <- function(coverage_ledger) {
  coverage_ledger[core_model_eligible %in% TRUE, feature]
}

s405m_order_age_levels <- function(values, cohort) {
  approved <- s405m_age_interval_levels(cohort)
  observed <- unique(trimws(as.character(values[!is.na(values)])))
  s405m_stop_if(!length(observed), "S405M_AGE_INTERVAL_NO_OBSERVED_LEVEL")
  s405m_stop_if(any(!observed %in% approved), "S405M_AGE_INTERVAL_UNAPPROVED_LEVEL")
  approved
}

s405m_encode_imputation_data <- function(x, cohort, eligible_fields) {
  analysis_fields <- unique(c(s405m_stay_fields(cohort), s405m_pair_fields()))
  eligible_fields <- intersect(analysis_fields, eligible_fields)
  s405m_stop_if(!length(eligible_fields), "S405M_NO_ELIGIBLE_COVARIATE")
  d <- as.data.frame(x[, c(
    "bias_spo2_minus_sao2", "paired_mean_saturation_c90", "sao2_anchor_icu_hours",
    eligible_fields
  ), with = FALSE])
  stay_key <- paste0(
    nchar(as.character(x$patient_key_internal)), ":", as.character(x$patient_key_internal), "|",
    nchar(as.character(x$stay_key_internal)), ":", as.character(x$stay_key_internal)
  )
  d$cluster_stay <- match(stay_key, unique(stay_key))
  d$cluster_patient <- as.integer(factor(x$patient_key_internal))
  metadata <- list(
    eligible_fields = eligible_fields,
    stay_fields = intersect(s405m_stay_fields(cohort), eligible_fields),
    pair_fields = intersect(s405m_pair_fields(), eligible_fields),
    sex_levels = c("female", "male"),
    age_levels = NULL,
    age_dummy_fields = character()
  )
  if ("sex_class" %in% eligible_fields) {
    sex <- s405m_normalize_sex(d$sex_class)
    d$sex_class <- ifelse(is.na(sex), NA_real_, as.numeric(sex == "male"))
  }
  if ("age_interval" %in% eligible_fields) {
    metadata$age_levels <- s405m_order_age_levels(d$age_interval, cohort)
    d$age_interval <- factor(
      trimws(as.character(d$age_interval)),
      levels = metadata$age_levels,
      ordered = FALSE
    )
    s405m_stop_if(anyNA(d$age_interval), "S405M_AGE_INTERVAL_MISSING_REQUIRES_SEPARATE_ORDINAL_IMPLEMENTATION")
    age_predictor_levels <- metadata$age_levels[-1L]
    metadata$age_dummy_fields <- sprintf(".age_interval_dummy_%02d", seq_along(age_predictor_levels))
    for (i in seq_along(age_predictor_levels)) {
      d[[metadata$age_dummy_fields[[i]]]] <- as.numeric(d$age_interval == age_predictor_levels[[i]])
    }
  }
  for (field in setdiff(eligible_fields, c("sex_class", "age_interval"))) {
    d[[field]] <- suppressWarnings(as.numeric(d[[field]]))
  }
  s405m_stop_if(anyNA(d$cluster_stay) || anyNA(d$cluster_patient), "S405M_CLUSTER_ID_MISSING")
  list(data = d, metadata = metadata)
}

s405m_imputation_spec <- function(encoded, cohort) {
  d <- encoded$data
  metadata <- encoded$metadata
  methods <- mice::make.method(d)
  methods[] <- ""
  targets <- metadata$eligible_fields[vapply(metadata$eligible_fields, function(field) anyNA(d[[field]]), logical(1L))]
  methods[targets] <- "ml.lmer"
  predictor_matrix <- mice::make.predictorMatrix(d)
  predictor_matrix[,] <- 0L
  predictor_pool <- c(
    "bias_spo2_minus_sao2", "paired_mean_saturation_c90", "sao2_anchor_icu_hours",
    metadata$eligible_fields
  )
  for (target in targets) {
    predictors <- setdiff(predictor_pool, target)
    predictor_matrix[target, predictors] <- 1L
  }
  predictor_matrix[, c("cluster_stay", "cluster_patient")] <- 0L
  predictor_matrix[c("cluster_stay", "cluster_patient"), ] <- 0L
  variables_levels <- setNames(rep("", ncol(d)), names(d))
  variables_levels[metadata$stay_fields] <- "cluster_stay"
  levels_id <- setNames(vector("list", length(targets)), targets)
  for (target in targets) {
    if (target %in% metadata$pair_fields) {
      levels_id[[target]] <- if (cohort == "eICU") c("cluster_stay", "cluster_patient") else c("cluster_stay")
    } else {
      levels_id[[target]] <- if (cohort == "eICU") c("cluster_patient") else character(0)
    }
  }
  model <- setNames(vector("list", length(targets)), targets)
  for (target in targets) model[[target]] <- if (target == "sex_class") "binary" else "pmm"
  list(
    method = methods, predictor_matrix = predictor_matrix,
    variables_levels = variables_levels, levels_id = levels_id,
    model = model, targets = targets
  )
}

s405m_complete_data <- function(mids, encoded, action) {
  d <- if (is.null(mids)) encoded$data else mice::complete(mids, action = action)
  metadata <- encoded$metadata
  if ("sex_class" %in% metadata$eligible_fields) {
    rounded <- round(as.numeric(d$sex_class))
    s405m_stop_if(any(!rounded %in% c(0, 1)), "S405M_IMPUTED_SEX_INVALID")
    d$sex_class <- factor(metadata$sex_levels[rounded + 1L], levels = metadata$sex_levels)
  }
  if ("age_interval" %in% metadata$eligible_fields) {
    label <- trimws(as.character(d$age_interval))
    s405m_stop_if(any(is.na(label) | !label %in% metadata$age_levels), "S405M_IMPUTED_AGE_INTERVAL_INVALID")
    d$age_interval <- factor(label, levels = metadata$age_levels, ordered = FALSE)
  }
  if (length(metadata$age_dummy_fields)) {
    d <- d[, setdiff(names(d), metadata$age_dummy_fields), drop = FALSE]
  }
  d
}

s405m_mids_logged_events <- function(mids, block, chain = NA_integer_) {
  if (is.null(mids) || is.null(mids$loggedEvents)) return(data.table::data.table())
  out <- data.table::as.data.table(mids$loggedEvents)
  out[, `:=`(imputation_block = block, external_chain = chain)]
  out
}

s405m_mids_trace <- function(mids, block, chain_offset = 0L) {
  if (is.null(mids) || is.null(mids$chainMean) || !length(mids$chainMean)) return(data.table::data.table())
  means <- mids$chainMean
  variances <- mids$chainVar
  dimensions <- dim(means)
  variables <- dimnames(means)[[1L]]
  rows <- list()
  for (variable_index in seq_len(dimensions[[1L]])) {
    for (iteration in seq_len(dimensions[[2L]])) {
      for (chain in seq_len(dimensions[[3L]])) {
        mean_value <- means[variable_index, iteration, chain]
        variance_value <- variances[variable_index, iteration, chain]
        if (!is.na(mean_value) || !is.na(variance_value)) {
          rows[[length(rows) + 1L]] <- data.table::data.table(
            imputation_block = block, feature = variables[[variable_index]],
            iteration = iteration, chain = chain_offset + chain,
            chain_mean = mean_value, chain_variance = variance_value
          )
        }
      }
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

s405m_trace_diagnostics <- function(trace) {
  if (!nrow(trace)) return(data.table::data.table())
  trace[, {
    chain_summary <- .SD[, .(
      chain_mean_average = mean(chain_mean, na.rm = TRUE),
      chain_mean_variance = stats::var(chain_mean, na.rm = TRUE),
      lag1_autocorrelation = if (.N >= 3L) suppressWarnings(stats::cor(chain_mean[-.N], chain_mean[-1L], use = "complete.obs")) else NA_real_
    ), by = chain]
    iterations <- uniqueN(iteration)
    within <- mean(chain_summary$chain_mean_variance, na.rm = TRUE)
    between <- iterations * stats::var(chain_summary$chain_mean_average, na.rm = TRUE)
    rhat <- if (is.finite(within) && within > 0 && is.finite(between)) {
      sqrt((((iterations - 1) / iterations) * within + between / iterations) / within)
    } else {
      NA_real_
    }
    .(
      chain_n = uniqueN(chain), iteration_n = iterations,
      pilot_rhat_screen = rhat,
      mean_absolute_lag1_autocorrelation = mean(abs(chain_summary$lag1_autocorrelation), na.rm = TRUE),
      final_chain_mean_min = min(chain_mean[iteration == max(iteration)], na.rm = TRUE),
      final_chain_mean_max = max(chain_mean[iteration == max(iteration)], na.rm = TRUE)
    )
  }, by = .(imputation_block, feature)]
}

s405m_predictor_matrix_ledger <- function(predictor_matrix, methods, block, chain = NA_integer_) {
  targets <- names(methods)[nzchar(methods)]
  if (!length(targets)) {
    return(data.table::data.table(
      imputation_block = character(), external_chain = integer(), target = character(),
      predictor = character(), matrix_code = integer(), predictor_role = character(),
      included = logical(), target_method = character()
    ))
  }
  data.table::rbindlist(lapply(targets, function(target) {
    predictor_names <- colnames(predictor_matrix)
    code <- as.integer(predictor_matrix[target, , drop = TRUE])
    data.table::data.table(
      imputation_block = block,
      external_chain = as.integer(chain),
      target = target,
      predictor = predictor_names,
      matrix_code = code,
      predictor_role = data.table::fcase(
        code == -2L, "cluster_identifier",
        code == 0L, "excluded",
        default = "fixed_auxiliary_predictor"
      ),
      included = code != 0L,
      target_method = methods[[target]],
      exclusion_reason = data.table::fcase(
        code != 0L, "",
        predictor_names == target, "target_itself",
        predictor_names %in% c("cluster_stay", "cluster_patient"), "protected_internal_identifier",
        default = "not_selected_by_approved_whitelist_or_incomplete_auxiliary"
      )
    )
  }))
}

s405m_stay_imputation_data <- function(encoded, cohort) {
  d <- data.table::as.data.table(encoded$data)
  stay_fields <- encoded$metadata$stay_fields
  pair_fields <- encoded$metadata$pair_fields
  stay_source_fields <- c(stay_fields, encoded$metadata$age_dummy_fields)
  base <- d[, lapply(.SD, function(v) {
    observed <- unique(v[!is.na(v)])
    s405m_stop_if(length(observed) > 1L, "S405M_ENCODED_STAY_VALUE_CONFLICT")
    if (!length(observed)) v[NA_integer_] else observed[[1L]]
  }), by = .(cluster_patient, cluster_stay), .SDcols = stay_source_fields]
  aggregate_fields <- c("bias_spo2_minus_sao2", "paired_mean_saturation_c90", "sao2_anchor_icu_hours", pair_fields)
  aggregates <- d[, c(
    lapply(.SD, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)),
    lapply(.SD, function(v) if (all(is.na(v))) NA_real_ else stats::sd(v, na.rm = TRUE))
  ), by = .(cluster_patient, cluster_stay), .SDcols = aggregate_fields]
  data.table::setnames(
    aggregates,
    c("cluster_patient", "cluster_stay", paste0(aggregate_fields, "_mean"), paste0(aggregate_fields, "_sd"))
  )
  merge(base, aggregates, by = c("cluster_patient", "cluster_stay"), all.x = TRUE, sort = FALSE)
}

s405m_run_stay_imputation <- function(encoded, cohort, m, maxit, seed, print_flag) {
  stay_data <- as.data.frame(s405m_stay_imputation_data(encoded, cohort))
  stay_fields <- encoded$metadata$stay_fields
  if ("sex_class" %in% stay_fields) {
    stay_data$sex_class <- factor(stay_data$sex_class, levels = c(0, 1))
  }
  targets <- stay_fields[vapply(stay_fields, function(field) anyNA(stay_data[[field]]), logical(1L))]
  if (!length(targets)) {
    if ("sex_class" %in% stay_fields) stay_data$sex_class <- as.numeric(as.character(stay_data$sex_class))
    empty_methods <- mice::make.method(stay_data)
    empty_methods[] <- ""
    empty_matrix <- mice::make.predictorMatrix(stay_data)
    empty_matrix[,] <- 0L
    return(list(
      mids = NULL,
      completed = replicate(m, stay_data, simplify = FALSE),
      logged_events = data.table::data.table(), trace = data.table::data.table(), targets = targets,
      methods = empty_methods, predictor_matrix = empty_matrix,
      predictor_ledger = s405m_predictor_matrix_ledger(empty_matrix, empty_methods, "stay")
    ))
  }
  methods <- mice::make.method(stay_data)
  methods[] <- ""
  predictor_matrix <- mice::make.predictorMatrix(stay_data)
  predictor_matrix[,] <- 0L
  predictors <- setdiff(names(stay_data), c("cluster_stay", "cluster_patient", "age_interval"))
  complete_predictors <- predictors[vapply(stay_data[predictors], function(v) !anyNA(v), logical(1L))]
  for (target in targets) predictor_matrix[target, setdiff(complete_predictors, target)] <- 1L
  predictor_matrix[, c("cluster_stay", "cluster_patient")] <- 0L
  predictor_matrix[c("cluster_stay", "cluster_patient"), ] <- 0L
  if (cohort == "eICU") {
    stay_data$.patient_mi <- match(stay_data$cluster_patient, unique(stay_data$cluster_patient))
    methods <- c(methods, .patient_mi = "")
    predictor_matrix <- mice::make.predictorMatrix(stay_data)
    predictor_matrix[,] <- 0L
    predictors <- setdiff(names(stay_data), c("cluster_stay", "cluster_patient", ".patient_mi", "age_interval"))
    complete_predictors <- predictors[vapply(stay_data[predictors], function(v) !anyNA(v), logical(1L))]
    for (target in targets) {
      predictor_matrix[target, setdiff(complete_predictors, target)] <- 1L
      predictor_matrix[target, ".patient_mi"] <- if (target == "sex_class") 0L else -2L
    }
    predictor_matrix[, c("cluster_stay", "cluster_patient")] <- 0L
    predictor_matrix[c("cluster_stay", "cluster_patient", ".patient_mi"), ] <- 0L
    methods[targets] <- ifelse(targets == "sex_class", "logreg", "2l.pmm.scaled")
    mids <- mice::mice(
      stay_data, m = m, maxit = maxit, method = methods,
      predictorMatrix = predictor_matrix, printFlag = print_flag,
      seed = seed, remove.collinear = FALSE, remove.constant = FALSE,
      donors = 5L
    )
  } else {
    methods[targets] <- ifelse(
      targets == "sex_class", "logreg",
      ifelse(targets == "age_interval", "polr", "pmm")
    )
    mids <- mice::mice(
      stay_data, m = m, maxit = maxit, method = methods,
      predictorMatrix = predictor_matrix, printFlag = print_flag,
      seed = seed, remove.collinear = FALSE, remove.constant = FALSE,
      donors = 5L
    )
  }
  completed <- lapply(seq_len(m), function(i) {
    out <- mice::complete(mids, action = i)
    if ("sex_class" %in% stay_fields) out$sex_class <- as.numeric(as.character(out$sex_class))
    out$.patient_mi <- NULL
    out
  })
  list(
    mids = mids,
    completed = completed,
    logged_events = s405m_mids_logged_events(mids, "stay"),
    trace = s405m_mids_trace(mids, "stay"), targets = targets,
    methods = methods, predictor_matrix = predictor_matrix,
    predictor_ledger = s405m_predictor_matrix_ledger(predictor_matrix, methods, "stay")
  )
}

s405m_run_pair_chain <- function(pair_data, cohort, pair_fields, maxit, seed, print_flag, chain) {
  targets <- pair_fields[vapply(pair_fields, function(field) anyNA(pair_data[[field]]), logical(1L))]
  if (!length(targets)) {
    empty_methods <- mice::make.method(pair_data)
    empty_methods[] <- ""
    empty_matrix <- mice::make.predictorMatrix(pair_data)
    empty_matrix[,] <- 0L
    return(list(
      mids = NULL, completed = pair_data,
      logged_events = data.table::data.table(), trace = data.table::data.table(), targets = targets,
      patient_auxiliary_n = 0L,
      patient_auxiliary_source_fields = character(),
      nonimputed_missing_predictors_excluded = character(),
      redundant_stay_auxiliaries_excluded = character(),
      methods = empty_methods, predictor_matrix = empty_matrix,
      predictor_ledger = s405m_predictor_matrix_ledger(empty_matrix, empty_methods, "pair", chain)
    ))
  }
  pair_data <- data.table::as.data.table(pair_data)
  pair_data[, .stay_mi := match(cluster_stay, unique(cluster_stay))]
  patient_auxiliary_n <- 0L
  patient_auxiliary_source_fields <- character()
  if (cohort == "eICU") {
    patient_auxiliary_source_fields <- intersect(
      c(
        "bias_spo2_minus_sao2", "paired_mean_saturation_c90", "sao2_anchor_icu_hours",
        pair_fields
      ),
      names(pair_data)
    )
    s405m_stop_if(
      any(!vapply(pair_data[, patient_auxiliary_source_fields, with = FALSE], is.numeric, logical(1L))),
      "S405M_PATIENT_AUXILIARY_NONNUMERIC_SOURCE"
    )
    patient_aux <- pair_data[, c(
      lapply(.SD, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)),
      list(.patient_stay_n = data.table::uniqueN(cluster_stay))
    ), by = cluster_patient, .SDcols = patient_auxiliary_source_fields]
    data.table::setnames(patient_aux, c("cluster_patient", paste0(".patient_aux_", patient_auxiliary_source_fields), ".patient_stay_n"))
    patient_auxiliary_n <- ncol(patient_aux) - 1L
    pair_data <- merge(pair_data, patient_aux, by = "cluster_patient", all.x = TRUE, sort = FALSE)
  }
  pair_data <- as.data.frame(pair_data)
  methods <- mice::make.method(pair_data)
  methods[] <- ""
  methods[targets] <- getOption("s405m.pair_imputation_method", "2l.pmm.scaled")
  predictor_matrix <- mice::make.predictorMatrix(pair_data)
  predictor_matrix[,] <- 0L
  approved_predictor_pool <- intersect(c(
    "bias_spo2_minus_sao2", "paired_mean_saturation_c90", "sao2_anchor_icu_hours",
    setdiff(s405m_stay_fields(cohort), "age_interval"),
    grep("^\\.age_interval_dummy_", names(pair_data), value = TRUE), pair_fields,
    paste0(".patient_aux_", patient_auxiliary_source_fields), ".patient_stay_n"
  ), names(pair_data))
  predictors <- approved_predictor_pool
  nonimputed_missing_predictors <- predictors[
    methods[predictors] == "" &
      vapply(pair_data[predictors], function(value) anyNA(value), logical(1L))
  ]
  redundant_stay_auxiliaries <- intersect(
    paste0(".patient_aux_", s405m_stay_fields(cohort)), predictors
  )
  usable_predictors <- setdiff(
    predictors,
    union(nonimputed_missing_predictors, redundant_stay_auxiliaries)
  )
  for (target in targets) {
    predictor_matrix[target, setdiff(usable_predictors, target)] <- 1L
    predictor_matrix[target, ".stay_mi"] <- -2L
  }
  predictor_matrix[, c("cluster_stay", "cluster_patient")] <- 0L
  predictor_matrix[c("cluster_stay", "cluster_patient", ".stay_mi"), ] <- 0L
  s405m_stop_if(any(vapply(nonimputed_missing_predictors, function(field) {
    any(predictor_matrix[targets, field, drop = FALSE] != 0L)
  }, logical(1L))), "S405M_NONIMPUTED_MISSING_PREDICTOR_RETAINED")
  s405m_stop_if(any(vapply(redundant_stay_auxiliaries, function(field) {
    any(predictor_matrix[targets, field, drop = FALSE] != 0L)
  }, logical(1L))), "S405M_REDUNDANT_STAY_AUXILIARY_RETAINED")
  mids <- mice::mice(
    pair_data, m = 1L, maxit = maxit, method = methods,
    predictorMatrix = predictor_matrix, printFlag = print_flag,
    seed = seed, remove.collinear = FALSE, remove.constant = FALSE,
    donors = 5L
  )
  completed <- mice::complete(mids, action = 1L)
  completed <- completed[, !grepl("^\\.patient_aux_|^\\.patient_stay_n$|^\\.stay_mi$", names(completed)), drop = FALSE]
  list(
    mids = mids, completed = completed,
    logged_events = s405m_mids_logged_events(mids, "pair", chain),
    trace = s405m_mids_trace(mids, "pair", chain_offset = chain - 1L),
    targets = targets, patient_auxiliary_n = patient_auxiliary_n,
    patient_auxiliary_source_fields = patient_auxiliary_source_fields,
    nonimputed_missing_predictors_excluded = nonimputed_missing_predictors,
    redundant_stay_auxiliaries_excluded = redundant_stay_auxiliaries,
    methods = methods, predictor_matrix = predictor_matrix,
    predictor_ledger = s405m_predictor_matrix_ledger(predictor_matrix, methods, "pair", chain)
  )
}

s405m_run_imputation <- function(x, cohort, eligible_fields, m = 3L, maxit = 5L,
                                 seed = 20260902L, print_flag = FALSE,
                                 execution_mode = c("pilot", "formal")) {
  execution_mode <- match.arg(execution_mode)
  if (execution_mode == "formal") {
    stop("S405M_FORMAL_CONVERGENCE_ADAPTER_NOT_APPROVED", call. = FALSE)
  }
  s405m_stop_if(m > 3L || maxit > 5L, "S405M_PILOT_LIMIT_EXCEEDED")
  prepared <- s405m_prepare_analytic_states(x, cohort)
  encoded <- s405m_encode_imputation_data(prepared$data, cohort, eligible_fields)
  stay_result <- s405m_run_stay_imputation(encoded, cohort, as.integer(m), as.integer(maxit), as.integer(seed), print_flag)
  raw <- data.table::as.data.table(encoded$data)
  pair_results <- vector("list", as.integer(m))
  completed <- vector("list", as.integer(m))
  for (i in seq_len(m)) {
    stay_completed <- data.table::as.data.table(stay_result$completed[[i]])
    stay_values <- stay_completed[, c("cluster_patient", "cluster_stay", encoded$metadata$stay_fields), with = FALSE]
    s405m_assert_unique(stay_values, c("cluster_patient", "cluster_stay"), "S405M_COMPLETED_STAY_COMPOSITE_KEY_DUPLICATE")
    pair_data <- merge(
      raw[, setdiff(names(raw), encoded$metadata$stay_fields), with = FALSE],
      stay_values, by = c("cluster_patient", "cluster_stay"), all.x = TRUE, sort = FALSE
    )
    s405m_stop_if(nrow(pair_data) != nrow(raw), "S405M_STAY_MAPPING_CARDINALITY")
    pair_data <- as.data.frame(pair_data[, names(encoded$data), with = FALSE])
    pair_results[[i]] <- s405m_run_pair_chain(
      pair_data, cohort, encoded$metadata$pair_fields,
      maxit = as.integer(maxit), seed = as.integer(seed + 10000L + i),
      print_flag = print_flag, chain = i
    )
    completed[[i]] <- s405m_complete_data(NULL, list(data = pair_results[[i]]$completed, metadata = encoded$metadata), 1L)
  }
  logged <- data.table::rbindlist(c(
    list(stay_result$logged_events),
    lapply(pair_results, function(x) x$logged_events)
  ), fill = TRUE)
  trace <- data.table::rbindlist(c(
    list(stay_result$trace),
    lapply(pair_results, function(x) x$trace)
  ), fill = TRUE)
  predictor_ledger <- data.table::rbindlist(c(
    list(stay_result$predictor_ledger),
    lapply(pair_results, function(x) x$predictor_ledger)
  ), fill = TRUE)
  sex_summary <- prepared$sex_ledger[, .(
    stay_n = .N,
    source_observed_n = sum(source_observed),
    patient_reconciled_n = sum(reconciled_from_patient),
    analytic_observed_before_imputation_n = sum(analytic_observed),
    residual_missing_n = sum(residual_missing)
  )]
  list(
    stay_mids = stay_result$mids, pair_mids = lapply(pair_results, function(x) x$mids),
    completed = completed, encoded = encoded,
    spec = list(
      stay_targets = stay_result$targets,
      pair_targets = unique(unlist(lapply(pair_results, function(x) x$targets))),
      pair_imputation_structure = if (cohort == "eICU") "pairs_within_stay_with_patient_shared_auxiliaries" else "pairs_within_patient_stay",
      patient_auxiliary_n = max(vapply(pair_results, function(x) x$patient_auxiliary_n, integer(1L))),
      patient_auxiliary_source_fields = unique(unlist(lapply(pair_results, function(x) x$patient_auxiliary_source_fields))),
      pair_nonimputed_missing_predictors_excluded = unique(unlist(lapply(pair_results, function(x) x$nonimputed_missing_predictors_excluded))),
      pair_redundant_stay_auxiliaries_excluded = unique(unlist(lapply(pair_results, function(x) x$redundant_stay_auxiliaries_excluded))),
      formal_convergence_ready = FALSE,
      convergence_interpretation = "pilot_interface_screen_only_not_formal_convergence_evidence"
    ),
    sex_reconciliation_summary = sex_summary,
    predictor_ledger = predictor_ledger,
    logged_events = logged, trace = trace,
    convergence_mean = if (!nrow(trace)) data.table::data.table() else s405m_trace_diagnostics(trace[!is.na(chain_mean)]),
    convergence_sd = if (!nrow(trace)) data.table::data.table() else s405m_trace_diagnostics(trace[!is.na(chain_variance), .(
      imputation_block, feature, iteration, chain,
      chain_mean = chain_variance, chain_variance = NA_real_
    )]),
    m = as.integer(m), maxit = as.integer(maxit)
  )
}

s405m_validate_imputation <- function(original_x, imputation_result, cohort) {
  encoded_original <- imputation_result$encoded$data
  completed <- imputation_result$completed
  fields <- imputation_result$encoded$metadata$eligible_fields
  rows <- list()
  for (field in fields) {
    original <- encoded_original[[field]]
    observed <- !is.na(original)
    observed_unchanged <- all(vapply(completed, function(d) {
      current <- d[[field]]
      if (field == "sex_class") current <- as.numeric(current == "male")
      if (field == "age_interval") current <- as.numeric(current)
      isTRUE(all.equal(as.numeric(current[observed]), as.numeric(original[observed]), tolerance = 0, check.attributes = FALSE))
    }, logical(1L)))
    completed_no_missing <- all(vapply(completed, function(d) !anyNA(d[[field]]), logical(1L)))
    rows[[length(rows) + 1L]] <- data.table::data.table(
      feature = field, observed_unchanged = observed_unchanged,
      completed_no_missing = completed_no_missing
    )
  }
  stay_fields <- intersect(s405m_stay_fields(cohort), fields)
  stay_constancy <- data.table::rbindlist(lapply(seq_along(completed), function(i) {
    d <- data.table::as.data.table(completed[[i]])
    if (!length(stay_fields)) return(data.table::data.table(imputation = i, feature = character(), stay_constant = logical()))
    data.table::rbindlist(lapply(stay_fields, function(field) {
      counts <- d[, .(value_n = data.table::uniqueN(get(field))), by = .(cluster_patient, cluster_stay)]
      data.table::data.table(imputation = i, feature = field, stay_constant = all(counts$value_n == 1L))
    }))
  }), fill = TRUE)
  patient_sex_constancy <- if (!"sex_class" %in% fields) {
    data.table::data.table(imputation = integer(), patient_constant = logical())
  } else {
    data.table::rbindlist(lapply(seq_along(completed), function(i) {
      d <- data.table::as.data.table(completed[[i]])
      counts <- d[, .(value_n = data.table::uniqueN(sex_class)), by = cluster_patient]
      data.table::data.table(imputation = i, patient_constant = all(counts$value_n == 1L))
    }))
  }
  list(
    variable = data.table::rbindlist(rows),
    stay_constancy = stay_constancy,
    patient_sex_constancy = patient_sex_constancy
  )
}

s405m_imputation_domain_spec <- function(field, metadata) {
  spec <- switch(field,
    age_value_per_10y = list(lower = 0, lower_inclusive = FALSE, upper = Inf, upper_inclusive = TRUE, integer_expected = FALSE),
    age_topcoded_indicator = list(lower = 0, lower_inclusive = TRUE, upper = 1, upper_inclusive = TRUE, integer_expected = TRUE),
    cci_value = list(lower = 0, lower_inclusive = TRUE, upper = Inf, upper_inclusive = TRUE, integer_expected = TRUE),
    sofa_total = list(lower = 0, lower_inclusive = TRUE, upper = 24, upper_inclusive = TRUE, integer_expected = TRUE),
    pao2_candidate_value = list(lower = 0, lower_inclusive = FALSE, upper = 1000, upper_inclusive = TRUE, integer_expected = FALSE),
    paco2_candidate_value = list(lower = 0, lower_inclusive = FALSE, upper = 1000, upper_inclusive = TRUE, integer_expected = FALSE),
    ph_candidate_value = list(lower = 6.3, lower_inclusive = TRUE, upper = 8.0, upper_inclusive = TRUE, integer_expected = FALSE),
    lactate_candidate_value = list(lower = 0, lower_inclusive = FALSE, upper = Inf, upper_inclusive = TRUE, integer_expected = FALSE),
    hb_candidate_value = list(lower = 0, lower_inclusive = FALSE, upper = 30, upper_inclusive = TRUE, integer_expected = FALSE),
    sex_class = list(lower = 0, lower_inclusive = TRUE, upper = 1, upper_inclusive = TRUE, integer_expected = TRUE),
    age_interval = list(lower = 1, lower_inclusive = TRUE, upper = length(metadata$age_levels), upper_inclusive = TRUE, integer_expected = TRUE),
    NULL
  )
  if (is.null(spec)) stop(paste0("S405M_IMPUTATION_DOMAIN_SPEC_MISSING_", field), call. = FALSE)
  spec
}

s405m_domain_numeric <- function(value, field, metadata) {
  if (field == "sex_class") {
    if (is.numeric(value)) return(as.numeric(value))
    label <- tolower(trimws(as.character(value)))
    return(ifelse(label == "female", 0, ifelse(label == "male", 1, NA_real_)))
  }
  if (field == "age_interval") {
    if (is.numeric(value)) return(as.numeric(value))
    return(match(trimws(as.character(value)), metadata$age_levels))
  }
  suppressWarnings(as.numeric(value))
}

s405m_imputed_domain_audit <- function(imputation_result) {
  original <- imputation_result$encoded$data
  metadata <- imputation_result$encoded$metadata
  fields <- metadata$eligible_fields
  data.table::rbindlist(lapply(seq_along(imputation_result$completed), function(imputation_index) {
    completed <- imputation_result$completed[[imputation_index]]
    data.table::rbindlist(lapply(fields, function(field) {
      missing_before <- is.na(original[[field]])
      value <- s405m_domain_numeric(completed[[field]], field, metadata)
      spec <- s405m_imputation_domain_spec(field, metadata)
      imputed <- value[missing_before]
      missing_after_n <- sum(is.na(imputed))
      nonfinite_n <- sum(!is.na(imputed) & !is.finite(imputed))
      finite <- !is.na(imputed) & is.finite(imputed)
      observed_support <- unique(s405m_domain_numeric(original[[field]][!missing_before], field, metadata))
      outside_support_n <- sum(finite & !imputed %in% observed_support)
      below_n <- sum(finite & if (spec$lower_inclusive) imputed < spec$lower else imputed <= spec$lower)
      above_n <- sum(finite & if (spec$upper_inclusive) imputed > spec$upper else imputed >= spec$upper)
      noninteger_n <- if (spec$integer_expected) sum(finite & abs(imputed - round(imputed)) > 1e-08) else 0L
      invalid_n <- missing_after_n + nonfinite_n + outside_support_n + below_n + above_n + noninteger_n
      data.table::data.table(
        imputation = as.integer(imputation_index), feature = field,
        imputed_cell_n = sum(missing_before), missing_after_n = missing_after_n,
        nonfinite_n = nonfinite_n, below_lower_bound_n = below_n,
        above_upper_bound_n = above_n, noninteger_support_warning_n = noninteger_n,
        outside_observed_support_n = outside_support_n,
        imputed_min = if (any(finite)) min(imputed[finite]) else NA_real_,
        imputed_max = if (any(finite)) max(imputed[finite]) else NA_real_,
        lower_bound = spec$lower, lower_inclusive = spec$lower_inclusive,
        upper_bound = spec$upper, upper_inclusive = spec$upper_inclusive,
        domain_pass = invalid_n == 0L
      )
    }))
  }))
}

s405m_imputation_missingness_manifest <- function(imputation_result) {
  original <- data.table::as.data.table(imputation_result$encoded$data)
  metadata <- imputation_result$encoded$metadata
  data.table::rbindlist(lapply(metadata$eligible_fields, function(field) {
    stay_level <- field %in% metadata$stay_fields
    if (stay_level) {
      unit <- original[, .(
        missing_state_n = data.table::uniqueN(is.na(get(field))),
        missing = all(is.na(get(field)))
      ), by = .(cluster_patient, cluster_stay)]
      s405m_stop_if(any(unit$missing_state_n != 1L), "S405M_STAY_MISSINGNESS_STATE_CONFLICT")
      unit_denominator_n <- nrow(unit)
      missing_unit_n <- sum(unit$missing)
    } else {
      unit_denominator_n <- nrow(original)
      missing_unit_n <- sum(is.na(original[[field]]))
    }
    data.table::data.table(
      feature = field,
      variable_level = if (stay_level) "study_icu_stay" else "pair",
      unit_denominator_n = unit_denominator_n,
      missing_unit_n = missing_unit_n,
      missing_pair_row_n = sum(is.na(original[[field]]))
    )
  }))
}

s405m_distribution_statistics <- function(value, prefix) {
  value <- value[!is.na(value) & is.finite(value)]
  names <- paste0(prefix, c(
    "_n", "_mean", "_sd", "_min", "_p01", "_p25", "_median", "_p75", "_p99", "_max"
  ))
  if (!length(value)) return(setNames(as.list(c(0, rep(NA_real_, 9L))), names))
  quantile_value <- stats::quantile(value, probs = c(0.01, 0.25, 0.5, 0.75, 0.99), names = FALSE, type = 7)
  setNames(as.list(c(
    length(value), mean(value), if (length(value) > 1L) stats::sd(value) else NA_real_,
    min(value), quantile_value, max(value)
  )), names)
}

s405m_observed_imputed_distribution <- function(imputation_result) {
  original <- data.table::as.data.table(imputation_result$encoded$data)
  metadata <- imputation_result$encoded$metadata
  data.table::rbindlist(lapply(metadata$eligible_fields, function(field) {
    stay_level <- field %in% metadata$stay_fields
    original_unit <- if (stay_level) {
      original[, .(value = get(field)[[1L]]), by = .(cluster_patient, cluster_stay)]
    } else {
      data.table::data.table(unit = seq_len(nrow(original)), value = original[[field]])
    }
    original_numeric <- s405m_domain_numeric(original_unit$value, field, metadata)
    observed <- original_numeric[!is.na(original_numeric)]
    imputed <- unlist(lapply(imputation_result$completed, function(completed) {
      completed <- data.table::as.data.table(completed)
      completed_unit <- if (stay_level) {
        completed[, .(value = get(field)[[1L]]), by = .(cluster_patient, cluster_stay)]
      } else {
        data.table::data.table(unit = seq_len(nrow(completed)), value = completed[[field]])
      }
      completed_numeric <- s405m_domain_numeric(completed_unit$value, field, metadata)
      completed_numeric[is.na(original_numeric)]
    }), use.names = FALSE)
    observed_stats <- s405m_distribution_statistics(observed, "observed")
    imputed_stats <- s405m_distribution_statistics(imputed, "imputed")
    observed_sd <- observed_stats$observed_sd
    standardized_mean_difference <- if (is.finite(observed_sd) && observed_sd > 0 && is.finite(imputed_stats$imputed_mean)) {
      (imputed_stats$imputed_mean - observed_stats$observed_mean) / observed_sd
    } else NA_real_
    sd_ratio <- if (is.finite(observed_sd) && observed_sd > 0 && is.finite(imputed_stats$imputed_sd)) {
      imputed_stats$imputed_sd / observed_sd
    } else NA_real_
    distribution_review_flag <-
      (is.finite(standardized_mean_difference) && abs(standardized_mean_difference) > 1) ||
      (is.finite(imputed_stats$imputed_median) && (
        imputed_stats$imputed_median < observed_stats$observed_p01 ||
          imputed_stats$imputed_median > observed_stats$observed_p99
      )) ||
      (is.finite(sd_ratio) && imputed_stats$imputed_n >= 20L && (sd_ratio < 0.1 || sd_ratio > 10))
    data.table::as.data.table(c(
      list(feature = field, variable_level = if (stay_level) "study_icu_stay" else "pair"),
      observed_stats, imputed_stats,
      list(
        standardized_mean_difference = standardized_mean_difference,
        imputed_to_observed_sd_ratio = sd_ratio,
        distribution_review_flag = distribution_review_flag
      )
    ))
  }), fill = TRUE)
}

s405m_fixed_formula <- function(cohort, predictor_fields) {
  random_term <- if (cohort == "eICU") "(1|cluster_patient/cluster_stay)" else "(1|cluster_patient)"
  stats::as.formula(paste0(
    "bias_spo2_minus_sao2 ~ ",
    paste(c("paired_mean_saturation_c90", predictor_fields), collapse = " + "),
    " + ", random_term
  ))
}

s405m_fit_lmms <- function(imputation_result, cohort, predictor_fields) {
  formula <- s405m_fixed_formula(cohort, predictor_fields)
  fits <- lapply(imputation_result$completed, function(d) {
    d$cluster_patient <- factor(d$cluster_patient)
    d$cluster_stay <- factor(d$cluster_stay)
    fit <- tryCatch(
      lme4::lmer(
        formula, data = d, REML = TRUE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L))
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      fit <- tryCatch(
        lme4::lmer(
          formula, data = d, REML = TRUE,
          control = lme4::lmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 100000L))
        ),
        error = function(e) e
      )
      s405m_stop_if(inherits(fit, "error"), "S405M_LMM_FIT_FAILURE")
      attr(fit, "optimizer_used") <- "nloptwrap"
    } else {
      attr(fit, "optimizer_used") <- "bobyqa"
    }
    fit
  })
  list(fits = fits, formula = formula, predictors = predictor_fields)
}

s405m_pool_estimates <- function(coefficients, variances, complete_df = Inf) {
  coefficients <- as.matrix(coefficients)
  variances <- as.matrix(variances)
  m <- ncol(coefficients)
  qbar <- rowMeans(coefficients)
  ubar <- rowMeans(variances)
  between <- apply(coefficients, 1L, stats::var)
  total <- ubar + (1 + 1 / m) * between
  relative_increase <- ifelse(ubar > 0, (1 + 1 / m) * between / ubar, 0)
  df_old <- ifelse(relative_increase > 0, (m - 1) * (1 + 1 / relative_increase)^2, Inf)
  lambda <- ifelse(total > 0, (1 + 1 / m) * between / total, 0)
  if (is.finite(complete_df)) {
    df_obs <- ((complete_df + 1) / (complete_df + 3)) * complete_df * (1 - lambda)
    df <- 1 / (1 / df_old + 1 / pmax(df_obs, 1e-08))
  } else {
    df <- df_old
  }
  se <- sqrt(total)
  critical <- ifelse(is.finite(df), stats::qt(0.975, df = df), stats::qnorm(0.975))
  statistic <- qbar / se
  p <- ifelse(is.finite(df), 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE), 2 * stats::pnorm(abs(statistic), lower.tail = FALSE))
  fmi <- (relative_increase + 2 / (df + 3)) / (relative_increase + 1)
  mcse <- sqrt(between / m)
  data.table::data.table(
    term = rownames(coefficients), estimate = qbar, std_error = se,
    ci_lower = qbar - critical * se, ci_upper = qbar + critical * se,
    statistic = statistic, df = df, p_value = p,
    within_variance = ubar, between_variance = between,
    fraction_missing_information = fmi,
    monte_carlo_error_estimate = mcse,
    monte_carlo_error_to_se = ifelse(se > 0, mcse / se, 0),
    imputations = m
  )
}

s405m_pool_lmms <- function(fit_result) {
  terms <- names(lme4::fixef(fit_result$fits[[1L]]))
  coefficients <- do.call(cbind, lapply(fit_result$fits, function(f) lme4::fixef(f)[terms]))
  variances <- do.call(cbind, lapply(fit_result$fits, function(f) diag(as.matrix(stats::vcov(f)))[terms]))
  patients <- length(unique(fit_result$fits[[1L]]@frame$cluster_patient))
  s405m_pool_estimates(coefficients, variances, complete_df = max(1, patients - length(terms)))
}

s405m_lmm_diagnostics <- function(fit_result) {
  fits <- fit_result$fits
  status <- data.table::rbindlist(lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    messages <- fit@optinfo$conv$lme4$messages
    vc <- data.table::as.data.table(as.data.frame(lme4::VarCorr(fit)))
    random <- vc[grp != "Residual" & var1 == "(Intercept)" & !is.na(vcov)]
    data.table::data.table(
      imputation = i, optimizer = attr(fit, "optimizer_used"),
      converged = is.null(messages), convergence_message = paste(messages, collapse = " | "),
      singular = lme4::isSingular(fit, tol = 1e-06),
      random_component_n = nrow(random),
      random_variance_min = if (nrow(random)) min(random$vcov) else NA_real_,
      random_variance_max = if (nrow(random)) max(random$vcov) else NA_real_
    )
  }))
  random_effects <- data.table::rbindlist(lapply(seq_along(fits), function(i) {
    vc <- data.table::as.data.table(as.data.frame(lme4::VarCorr(fits[[i]])))
    vc[, imputation := i]
    vc
  }), fill = TRUE)
  list(status = status, random_effects = random_effects)
}

s405m_fit_patient_balanced <- function(imputation_result, predictor_fields) {
  fixed_formula <- stats::as.formula(paste0(
    "bias_spo2_minus_sao2 ~ ",
    paste(c("paired_mean_saturation_c90", predictor_fields), collapse = " + ")
  ))
  fits <- lapply(imputation_result$completed, function(d) {
    patient_n <- table(d$cluster_patient)
    d$patient_balance_weight <- 1 / as.numeric(patient_n[as.character(d$cluster_patient)])
    fit <- stats::lm(fixed_formula, data = d, weights = patient_balance_weight)
    covariance <- sandwich::vcovCL(fit, cluster = d$cluster_patient, type = "HC1")
    list(fit = fit, covariance = covariance)
  })
  terms <- names(stats::coef(fits[[1L]]$fit))
  coefficients <- do.call(cbind, lapply(fits, function(x) stats::coef(x$fit)[terms]))
  variances <- do.call(cbind, lapply(fits, function(x) diag(x$covariance)[terms]))
  patient_n <- length(unique(imputation_result$completed[[1L]]$cluster_patient))
  pooled <- s405m_pool_estimates(coefficients, variances, complete_df = max(1, patient_n - length(terms)))
  list(pooled = pooled, patient_n = patient_n, fits = fits)
}

s405m_compare_primary_balanced <- function(primary, balanced) {
  out <- merge(
    primary[, .(term, primary_estimate = estimate, primary_ci_lower = ci_lower, primary_ci_upper = ci_upper)],
    balanced[, .(term, balanced_estimate = estimate, balanced_ci_lower = ci_lower, balanced_ci_upper = ci_upper)],
    by = "term", all = TRUE
  )
  out[, same_direction := sign(primary_estimate) == sign(balanced_estimate)]
  out[, estimate_difference := balanced_estimate - primary_estimate]
  out
}
