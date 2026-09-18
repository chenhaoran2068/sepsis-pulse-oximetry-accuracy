r211b_pair_key_columns <- function() {
  c("patient_key_internal", "stay_key_internal", "sao2_stable_event_ordinal")
}

r211b_expected_metric_units <- function() {
  c(pao2 = "mmHg", paco2 = "mmHg", ph = "pH", lactate = "mmol/L", hb = "g/dL")
}

r211b_require_columns <- function(x, required, code) {
  if (!all(required %in% names(x))) stop(code, call. = FALSE)
  invisible(TRUE)
}

r211b_safe_error_code <- function(condition) {
  message <- conditionMessage(condition)
  if (grepl("^R211B_[A-Z0-9_]+$", message)) return(message)
  "R211B_UNCLASSIFIED_RUNTIME_ERROR"
}

r211b_validate_candidate_root <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("R211B_CANDIDATE_ROOT_EXISTS", call. = FALSE)
  }
  invisible(TRUE)
}

r211b_validate_qa_root <- function(path) {
  if (!dir.exists(path)) return(invisible(TRUE))
  root_entries <- list.files(path, all.files = TRUE, no.. = TRUE)
  if (length(setdiff(root_entries, "runtime")) > 0L) stop("R211B_QA_ROOT_EXISTS", call. = FALSE)
  runtime_dir <- file.path(path, "runtime")
  if (dir.exists(runtime_dir)) {
    runtime_entries <- list.files(runtime_dir, all.files = TRUE, no.. = TRUE)
    allowed <- c("runner_stdout.log", "runner_stderr.log")
    if (length(setdiff(runtime_entries, allowed)) > 0L) stop("R211B_QA_ROOT_EXISTS", call. = FALSE)
  }
  invisible(TRUE)
}

r211b_key_string <- function(x, keys = r211b_pair_key_columns()) {
  do.call(paste, c(x[, keys, drop = FALSE], sep = "\r"))
}

r211b_filter_m60 <- function(x, code) {
  r211b_require_columns(x, c(r211b_pair_key_columns(), "in_m60"), paste0(code, "_M60_INTERFACE_MISSING"))
  keep <- !is.na(x$in_m60) & as.logical(x$in_m60)
  y <- x[keep, , drop = FALSE]
  if (anyDuplicated(y[, r211b_pair_key_columns(), drop = FALSE])) stop(paste0(code, "_M60_DUPLICATE_KEY"), call. = FALSE)
  y
}

r211b_normalize_r07_metric_schema <- function(x) {
  source_prefix <- c("PaO2", "PaCO2", "pH", "Lactate", "Hb")
  target_prefix <- c("pao2", "paco2", "ph", "lactate", "hb")
  suffixes <- c("candidate_value", "canonical_unit", "linkage_class", "value_state")
  source_names <- as.vector(outer(source_prefix, suffixes, paste, sep = "_"))
  target_names <- as.vector(outer(target_prefix, suffixes, paste, sep = "_"))
  present <- source_names %in% names(x)
  if (any(present & target_names %in% names(x))) stop("R211B_R07_METRIC_AMBIGUOUS_SCHEMA", call. = FALSE)
  if (any(present)) names(x)[match(source_names[present], names(x))] <- target_names[present]
  x
}

r211b_align_exact <- function(pairs, sidecar, code) {
  keys <- r211b_pair_key_columns()
  if (anyDuplicated(pairs[, keys, drop = FALSE])) stop(paste0(code, "_PAIR_DUPLICATE_KEY"), call. = FALSE)
  if (anyDuplicated(sidecar[, keys, drop = FALSE])) stop(paste0(code, "_SIDECAR_DUPLICATE_KEY"), call. = FALSE)
  pair_key <- r211b_key_string(pairs)
  sidecar_key <- r211b_key_string(sidecar)
  if (any(!pair_key %in% sidecar_key) || any(!sidecar_key %in% pair_key)) stop(paste0(code, "_KEY_COVERAGE_FAILURE"), call. = FALSE)
  sidecar[match(pair_key, sidecar_key), , drop = FALSE]
}

r211b_align_sofa <- function(pairs, sofa) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("R211B_DATATABLE_REQUIRED", call. = FALSE)
  required_pair <- c("pair_key_internal", "stay_key_internal", "sao2_relative_icu_minutes")
  required_sofa <- c("stay_key_internal", "hour_end_min", "sofa_total", "sofa_status")
  r211b_require_columns(pairs, required_pair, "R211B_PAIR_SOFA_INTERFACE_MISSING")
  r211b_require_columns(sofa, required_sofa, "R211B_HOURLY_SOFA_INTERFACE_MISSING")
  p <- data.table::as.data.table(pairs[, required_pair, drop = FALSE])
  s <- data.table::as.data.table(sofa[, required_sofa, drop = FALSE])
  if (anyDuplicated(s[, c("stay_key_internal", "hour_end_min"), with = FALSE])) stop("R211B_HOURLY_SOFA_DUPLICATE_STAY_HOUR", call. = FALSE)
  data.table::setkey(s, stay_key_internal, hour_end_min)
  aligned <- s[p,
    on = .(stay_key_internal, hour_end_min = sao2_relative_icu_minutes),
    roll = Inf,
    .(pair_key_internal = i.pair_key_internal,
      sofa_total = x.sofa_total,
      sofa_status = x.sofa_status,
      sofa_delta_minutes = i.sao2_relative_icu_minutes - x.hour_end_min)
  ]
  unavailable <- is.na(aligned$sofa_delta_minutes) | aligned$sofa_delta_minutes < 0 | aligned$sofa_delta_minutes > 120
  aligned$sofa_total[unavailable] <- NA_real_
  aligned$sofa_status[unavailable] <- "unavailable_under_anchor_rule"
  out <- as.data.frame(aligned)
  out
}

r211b_metric_columns <- function() {
  metrics <- names(r211b_expected_metric_units())
  unlist(lapply(metrics, function(metric) {
    c(paste0(metric, "_candidate_value"), paste0(metric, "_canonical_unit"), paste0(metric, "_linkage_class"), paste0(metric, "_value_state"))
  }), use.names = FALSE)
}

r211b_validate_metric_values <- function(metric) {
  required <- c(r211b_pair_key_columns(), "in_m60", r211b_metric_columns())
  r211b_require_columns(metric, required, "R211B_R07_METRIC_INTERFACE_MISSING")
  expected_units <- r211b_expected_metric_units()
  permitted_states <- c(
    "eligible_pending_candidate_build", "unavailable", "ambiguous",
    "numeric_parse_hold", "source_scale_hold",
    "nonpositive_invalid", "outside_broad_technical_range"
  )
  held_states <- setdiff(permitted_states, "eligible_pending_candidate_build")
  for (name in names(expected_units)) {
    value <- metric[[paste0(name, "_candidate_value")]]
    unit <- as.character(metric[[paste0(name, "_canonical_unit")]])
    state <- as.character(metric[[paste0(name, "_value_state")]])
    if (any(!is.na(state) & !state %in% permitted_states)) stop("R211B_R07_METRIC_STATE_UNEXPECTED", call. = FALSE)
    if (any(!is.na(value) & unit != expected_units[[name]])) stop("R211B_R07_METRIC_UNIT_INVALID", call. = FALSE)
    if (any(state == "eligible_pending_candidate_build" & is.na(value))) stop("R211B_R07_METRIC_ELIGIBLE_VALUE_ABSENT", call. = FALSE)
    if (any(state %in% held_states & !is.na(value))) stop("R211B_R07_METRIC_HELD_VALUE_PRESENT", call. = FALSE)
  }
  invisible(TRUE)
}

r211b_prepare_clinical <- function(pairs, clinical, cohort, column, permitted) {
  field <- paste0(column, "_state")
  scope_field <- paste0(column, "_scope")
  n <- nrow(pairs)
  if (!permitted) {
    return(data.frame(
      value = rep(NA_character_, n),
      scope = rep("excluded_by_gate_r08", n),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (is.null(clinical)) stop("R211B_R08_REQUIRED_INPUT_ABSENT", call. = FALSE)
  filtered <- r211b_filter_m60(clinical, paste0("R211B_", cohort, "_R08"))
  aligned <- r211b_align_exact(pairs, filtered, paste0("R211B_", cohort, "_R08"))
  r211b_require_columns(aligned, column, "R211B_R08_CLINICAL_COLUMN_MISSING")
  value <- tolower(trimws(as.character(aligned[[column]])))
  if (any(is.na(value) | !value %in% c("yes", "unknown"))) stop("R211B_R08_CLINICAL_STATE_INVALID", call. = FALSE)
  data.frame(value = value, scope = rep("permitted_by_gate_r08", n), stringsAsFactors = FALSE, check.names = FALSE)
}

r211b_feature_role_manifest <- function() {
  data.frame(
    field = c(
      "cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid",
      "bias_spo2_minus_sao2", "spo2_percent", "sao2_percent",
      "age_years_numeric", "age_interval", "age_model_representation", "sex_class", "cci_value", "sofa_total", "sofa_status",
      paste0(names(r211b_expected_metric_units()), "_candidate_value"),
      unlist(lapply(names(r211b_expected_metric_units()), function(x) c(paste0(x, "_canonical_unit"), paste0(x, "_linkage_class"), paste0(x, "_value_state")))),
      "concurrent_vasoactive_medication_use_state", "concurrent_vasoactive_medication_use_scope",
      "current_invasive_mechanical_ventilation_state", "current_invasive_mechanical_ventilation_scope"
    ),
    role = c(
      "protected", "protected", "protected", "protected",
      "outcome", "candidate_predictor", "candidate_predictor",
      "candidate_predictor", "candidate_predictor", "metadata", "candidate_predictor", "candidate_predictor", "candidate_predictor", "metadata",
      rep("candidate_predictor", 5L), rep("metadata", 15L),
      "candidate_predictor", "metadata", "candidate_predictor", "metadata"
    ),
    model_rule = c(
      "never_predictor", "never_predictor", "never_predictor", "never_predictor",
      "outcome_only", "observed_never_impute", "observed_never_impute",
      "use_only_when_source_numeric", "mutually_exclusive_with_age_interval", "source_faithful", "candidate_cohort_local", "candidate_cohort_local", "candidate_cohort_local", "never_predictor",
      rep("candidate_cohort_local", 5L), rep("never_predictor", 15L),
      "only_if_scope_permitted_and_state_yes_unknown", "never_predictor", "only_if_scope_permitted_and_state_yes_unknown", "never_predictor"
    ),
    stringsAsFactors = FALSE
  )
}

r211b_build_candidate <- function(cohort, pairs, core, sofa, metric, clinical = NULL,
                                  allow_vaso = FALSE, allow_imv = FALSE) {
  pair_required <- c("pair_key_internal", "cohort", r211b_pair_key_columns(), "spo2_saturation_percent", "sao2_saturation_percent", "sao2_relative_icu_minutes", "pair_window_minutes", "contains_possible_transient_artifact")
  core_required <- c("patient_key_internal", "stay_key_internal", "age_years_numeric", "age_interval", "age_model_representation", "sex_class", "cci_value")
  r211b_require_columns(pairs, pair_required, "R211B_R06_PAIR_INTERFACE_MISSING")
  r211b_require_columns(core, core_required, "R211B_R06_CORE_INTERFACE_MISSING")
  if (anyDuplicated(pairs$pair_key_internal) || anyNA(pairs$pair_key_internal)) stop("R211B_PAIR_UID_INTEGRITY_FAILURE", call. = FALSE)
  if (any(pairs$cohort != cohort) || any(pairs$pair_window_minutes != 60L) || any(pairs$contains_possible_transient_artifact)) stop("R211B_R06_M60_SCOPE_INVALID", call. = FALSE)
  if (any(!is.finite(pairs$spo2_saturation_percent)) || any(!is.finite(pairs$sao2_saturation_percent)) || any(pairs$spo2_saturation_percent < 70 | pairs$spo2_saturation_percent > 100) || any(pairs$sao2_saturation_percent < 70 | pairs$sao2_saturation_percent > 100)) stop("R211B_ACCEPTED_PAIR_SATURATION_INVALID", call. = FALSE)

  core_keys <- c("patient_key_internal", "stay_key_internal")
  if (anyDuplicated(core[, core_keys, drop = FALSE])) stop("R211B_CORE_DUPLICATE_STAY_KEY", call. = FALSE)
  pair_core_key <- r211b_key_string(pairs, core_keys)
  core_key <- r211b_key_string(core, core_keys)
  if (any(!pair_core_key %in% core_key)) stop("R211B_R06_PAIR_TO_CORE_UNMATCHED", call. = FALSE)
  core_aligned <- core[match(pair_core_key, core_key), , drop = FALSE]

  metric <- r211b_normalize_r07_metric_schema(metric)
  metric <- r211b_filter_m60(metric, paste0("R211B_", cohort, "_R07"))
  r211b_validate_metric_values(metric)
  metric_aligned <- r211b_align_exact(pairs, metric, paste0("R211B_", cohort, "_R07"))

  sofa_aligned <- r211b_align_sofa(pairs, sofa)
  sofa_aligned <- sofa_aligned[match(pairs$pair_key_internal, sofa_aligned$pair_key_internal), , drop = FALSE]
  if (anyDuplicated(sofa_aligned$pair_key_internal) || !identical(as.character(sofa_aligned$pair_key_internal), as.character(pairs$pair_key_internal))) stop("R211B_SOFA_ALIGNMENT_CARDINALITY_FAILURE", call. = FALSE)

  vaso <- r211b_prepare_clinical(pairs, clinical, cohort, "concurrent_vasoactive_medication_use", allow_vaso)
  imv <- r211b_prepare_clinical(pairs, clinical, cohort, "current_invasive_mechanical_ventilation", allow_imv)

  candidate <- data.frame(
    cohort = rep(cohort, nrow(pairs)),
    analysis_patient_key = pairs$patient_key_internal,
    analysis_stay_key = pairs$stay_key_internal,
    pair_uid = pairs$pair_key_internal,
    bias_spo2_minus_sao2 = pairs$spo2_saturation_percent - pairs$sao2_saturation_percent,
    spo2_percent = pairs$spo2_saturation_percent,
    sao2_percent = pairs$sao2_saturation_percent,
    age_years_numeric = core_aligned$age_years_numeric,
    age_interval = core_aligned$age_interval,
    age_model_representation = core_aligned$age_model_representation,
    sex_class = core_aligned$sex_class,
    cci_value = core_aligned$cci_value,
    sofa_total = sofa_aligned$sofa_total,
    sofa_status = sofa_aligned$sofa_status,
    concurrent_vasoactive_medication_use_state = vaso$value,
    concurrent_vasoactive_medication_use_scope = vaso$scope,
    current_invasive_mechanical_ventilation_state = imv$value,
    current_invasive_mechanical_ventilation_scope = imv$scope,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (name in names(r211b_expected_metric_units())) {
    candidate[[paste0(name, "_candidate_value")]] <- metric_aligned[[paste0(name, "_candidate_value")]]
    candidate[[paste0(name, "_canonical_unit")]] <- metric_aligned[[paste0(name, "_canonical_unit")]]
    candidate[[paste0(name, "_linkage_class")]] <- metric_aligned[[paste0(name, "_linkage_class")]]
    candidate[[paste0(name, "_value_state")]] <- metric_aligned[[paste0(name, "_value_state")]]
  }
  r211b_validate_candidate(candidate, cohort)
  candidate
}

r211b_validate_candidate <- function(candidate, cohort = NULL) {
  required <- r211b_feature_role_manifest()$field
  r211b_require_columns(candidate, required, "R211B_CANDIDATE_SCHEMA_INVALID")
  if (anyNA(candidate$analysis_patient_key) || anyNA(candidate$analysis_stay_key) || anyNA(candidate$pair_uid) || anyDuplicated(candidate$pair_uid)) stop("R211B_CANDIDATE_KEY_INTEGRITY_FAILURE", call. = FALSE)
  if (!is.null(cohort) && any(candidate$cohort != cohort)) stop("R211B_COHORT_SCOPE_INVALID", call. = FALSE)
  if (any(!is.finite(candidate$spo2_percent)) || any(!is.finite(candidate$sao2_percent)) || any(candidate$spo2_percent < 70 | candidate$spo2_percent > 100) || any(candidate$sao2_percent < 70 | candidate$sao2_percent > 100)) stop("R211B_CANDIDATE_SATURATION_INVALID", call. = FALSE)
  if (!isTRUE(all.equal(candidate$bias_spo2_minus_sao2, candidate$spo2_percent - candidate$sao2_percent, check.attributes = FALSE))) stop("R211B_BIAS_DIRECTION_INVALID", call. = FALSE)
  prohibited <- grep("relative_icu|lag|direction|artifact|source_label|source_code", names(candidate), value = TRUE)
  if (length(prohibited) > 0L) stop("R211B_PROTECTED_TIMING_OR_QC_FIELD_PRESENT", call. = FALSE)
  for (name in names(r211b_expected_metric_units())) {
    value <- candidate[[paste0(name, "_candidate_value")]]
    unit <- as.character(candidate[[paste0(name, "_canonical_unit")]])
    state <- as.character(candidate[[paste0(name, "_value_state")]])
    if (any(!is.na(value) & unit != r211b_expected_metric_units()[[name]])) stop("R211B_CANDIDATE_METRIC_UNIT_INVALID", call. = FALSE)
    if (any(state == "eligible_pending_candidate_build" & is.na(value))) stop("R211B_CANDIDATE_METRIC_STATE_VALUE_MISMATCH", call. = FALSE)
    if (any(state != "eligible_pending_candidate_build" & !is.na(value))) stop("R211B_CANDIDATE_METRIC_STATE_VALUE_MISMATCH", call. = FALSE)
    if (name == "ph" && any(!is.na(value) & (value < 6.3 | value > 8.0))) stop("R211B_CANDIDATE_PH_RANGE_INVALID", call. = FALSE)
  }
  for (field in c("concurrent_vasoactive_medication_use", "current_invasive_mechanical_ventilation")) {
    state <- candidate[[paste0(field, "_state")]]
    scope <- candidate[[paste0(field, "_scope")]]
    if (any(scope == "permitted_by_gate_r08" & (is.na(state) | !state %in% c("yes", "unknown")))) stop("R211B_CLINICAL_STATE_SCOPE_MISMATCH", call. = FALSE)
    if (any(scope == "excluded_by_gate_r08" & !is.na(state))) stop("R211B_CLINICAL_STATE_EXCLUDED_VALUE_PRESENT", call. = FALSE)
  }
  invisible(TRUE)
}

r211b_feature_completeness <- function(candidate) {
  manifest <- r211b_feature_role_manifest()
  features <- manifest$field[manifest$role == "candidate_predictor" & !grepl("_state$", manifest$field)]
  do.call(rbind, lapply(features, function(feature) {
    observed <- if (feature == "sex_class") {
      tolower(trimws(as.character(candidate[[feature]]))) %in% c("female", "male")
    } else {
      !is.na(candidate[[feature]])
    }
    data.frame(
      cohort = unique(candidate$cohort),
      feature = feature,
      eligible_pair_n = nrow(candidate),
      observed_pair_n = sum(observed),
      unavailable_or_unknown_pair_n = sum(!observed),
      stringsAsFactors = FALSE
    )
  }))
}

r211b_clinical_state_summary <- function(candidate) {
  do.call(rbind, lapply(c("concurrent_vasoactive_medication_use", "current_invasive_mechanical_ventilation"), function(feature) {
    state <- as.character(candidate[[paste0(feature, "_state")]])
    scope <- as.character(candidate[[paste0(feature, "_scope")]])
    permitted <- scope == "permitted_by_gate_r08"
    excluded <- scope == "excluded_by_gate_r08"
    if (any(!permitted & !excluded) || any(permitted & (is.na(state) | !state %in% c("yes", "unknown"))) || any(excluded & !is.na(state))) stop("R211B_CLINICAL_STATE_SUMMARY_INVALID", call. = FALSE)
    data.frame(
      cohort = unique(candidate$cohort),
      feature = feature,
      eligible_pair_n = nrow(candidate),
      permitted_pair_n = sum(permitted),
      yes_pair_n = sum(permitted & state == "yes"),
      unknown_pair_n = sum(permitted & state == "unknown"),
      excluded_by_gate_pair_n = sum(excluded),
      stringsAsFactors = FALSE
    )
  }))
}
