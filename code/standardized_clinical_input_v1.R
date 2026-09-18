ci1_stop <- function(code) stop(code, call. = FALSE)

ci1_stay_columns <- function() c(
  "patient_key_internal", "stay_key_internal", "age_years_numeric",
  "age_interval", "age_model_representation", "sex_class", "cci_value"
)

ci1_pair_keys <- function() c(
  "patient_key_internal", "stay_key_internal",
  "spo2_stable_event_ordinal", "sao2_stable_event_ordinal"
)

ci1_metric_names <- function() c("pao2", "paco2", "ph", "lactate", "hb")

ci1_metric_columns <- function() unlist(lapply(ci1_metric_names(), function(x) paste0(
  x, c("_candidate_value", "_canonical_unit", "_linkage_class", "_value_state")
)), use.names = FALSE)

ci1_treatment_names <- function() c(
  "concurrent_vasoactive_medication_use", "current_invasive_mechanical_ventilation"
)

ci1_pair_columns <- function() c(
  ci1_pair_keys(), "pair_key_internal", "sofa_total", "sofa_status",
  ci1_metric_columns(),
  unlist(lapply(ci1_treatment_names(), function(x) paste0(x, c("_state", "_scope"))), use.names = FALSE)
)

ci1_read_table <- function(path, expected) {
  if (length(path) != 1L || is.na(path) || !file.exists(path) || dir.exists(path)) ci1_stop("CI1_INPUT_FILE_ABSENT")
  first_line <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (!length(first_line)) ci1_stop("CI1_INPUT_EMPTY")
  header <- strsplit(first_line, "\t", fixed = TRUE)[[1L]]
  if (!identical(header, expected)) ci1_stop("CI1_INPUT_HEADER_MISMATCH")
  x <- tryCatch(utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, colClasses = "character", na.strings = "",
    fileEncoding = "UTF-8", stringsAsFactors = FALSE
  ), error = function(e) ci1_stop("CI1_INPUT_PARSE_ERROR"))
  if (!nrow(x)) ci1_stop("CI1_INPUT_EMPTY")
  if (!identical(names(x), expected)) ci1_stop("CI1_INPUT_HEADER_MISMATCH")
  x
}

ci1_number <- function(x, allow_missing = TRUE, integer = FALSE, code = "CI1_NUMBER_INVALID") {
  value <- suppressWarnings(as.numeric(x))
  bad <- (!is.na(x) & (!is.finite(value) | is.na(value))) |
    (!allow_missing & is.na(value)) |
    (integer & !is.na(value) & abs(value - round(value)) > 1e-8)
  if (any(bad)) ci1_stop(code)
  value
}

ci1_key <- function(x, fields) do.call(paste, c(x[fields], sep = "\r"))

ci1_read_stays <- function(path, eligible_stays) {
  x <- ci1_read_table(path, ci1_stay_columns())
  if (anyNA(x$patient_key_internal) || anyNA(x$stay_key_internal) ||
      any(!nzchar(x$patient_key_internal)) || any(!nzchar(x$stay_key_internal)) ||
      anyDuplicated(x$stay_key_internal)) ci1_stop("CI1_STAY_KEY_INVALID")
  if (!is.data.frame(eligible_stays) ||
      !all(c("patient_key_internal", "stay_key_internal") %in% names(eligible_stays)) ||
      anyDuplicated(eligible_stays$stay_key_internal)) ci1_stop("CI1_ELIGIBLE_STAY_INTERFACE_INVALID")
  if (!setequal(ci1_key(x, c("patient_key_internal", "stay_key_internal")),
                ci1_key(eligible_stays, c("patient_key_internal", "stay_key_internal"))))
    ci1_stop("CI1_STAY_COVERAGE_MISMATCH")
  x$age_years_numeric <- ci1_number(x$age_years_numeric, code = "CI1_AGE_INVALID")
  x$cci_value <- ci1_number(x$cci_value, integer = TRUE, code = "CI1_CCI_INVALID")
  if (any(!is.na(x$age_years_numeric) & x$age_years_numeric < 18) ||
      any(!is.na(x$cci_value) & x$cci_value < 0)) ci1_stop("CI1_STAY_VALUE_INVALID")
  if (any(!is.na(x$age_years_numeric) & !is.na(x$age_interval))) ci1_stop("CI1_AGE_REPRESENTATION_CONFLICT")
  if (any(!is.na(x$sex_class) & !x$sex_class %in% c("female", "male", "unknown_or_missing")))
    ci1_stop("CI1_SEX_INVALID")
  if (any(!is.na(x$age_model_representation) &
          !x$age_model_representation %in% c("numeric", "source_interval")))
    ci1_stop("CI1_AGE_REPRESENTATION_INVALID")
  x
}

ci1_read_pairs <- function(path, matched_pairs, stays, cohort, window_minutes) {
  x <- ci1_read_table(path, ci1_pair_columns())
  if (!is.data.frame(matched_pairs) || !all(c(ci1_pair_keys(), "cohort", "pair_window_minutes") %in% names(matched_pairs)))
    ci1_stop("CI1_MATCHED_PAIR_INTERFACE_INVALID")
  if (length(cohort) != 1L || length(window_minutes) != 1L ||
      any(matched_pairs$cohort != cohort) || any(matched_pairs$pair_window_minutes != window_minutes))
    ci1_stop("CI1_MATCHED_PAIR_SCOPE_INVALID")
  if (anyNA(x[c(ci1_pair_keys(), "pair_key_internal")]) ||
      any(vapply(x[c(ci1_pair_keys(), "pair_key_internal")], function(v) any(!nzchar(v)), logical(1L))) ||
      anyDuplicated(x$pair_key_internal)) ci1_stop("CI1_PAIR_KEY_INVALID")
  for (field in c("spo2_stable_event_ordinal", "sao2_stable_event_ordinal")) {
    x[[field]] <- ci1_number(x[[field]], allow_missing = FALSE, integer = TRUE, code = "CI1_PAIR_ORDINAL_INVALID")
    if (any(x[[field]] < 1)) ci1_stop("CI1_PAIR_ORDINAL_INVALID")
  }
  if (anyDuplicated(ci1_key(x, ci1_pair_keys())) || anyDuplicated(ci1_key(matched_pairs, ci1_pair_keys())) ||
      !setequal(ci1_key(x, ci1_pair_keys()), ci1_key(matched_pairs, ci1_pair_keys())))
    ci1_stop("CI1_PAIR_COVERAGE_MISMATCH")
  if (!all(ci1_key(x, c("patient_key_internal", "stay_key_internal")) %in%
           ci1_key(stays, c("patient_key_internal", "stay_key_internal"))))
    ci1_stop("CI1_PAIR_STAY_MISMATCH")
  x$sofa_total <- ci1_number(x$sofa_total, integer = TRUE, code = "CI1_SOFA_INVALID")
  if (any(!is.na(x$sofa_total) & (x$sofa_total < 0 | x$sofa_total > 24))) ci1_stop("CI1_SOFA_INVALID")
  if (any(is.na(x$sofa_status) | !nzchar(x$sofa_status))) ci1_stop("CI1_SOFA_STATUS_INVALID")
  if (any((x$sofa_status == "available") != !is.na(x$sofa_total))) ci1_stop("CI1_SOFA_STATE_VALUE_MISMATCH")
  units <- c(pao2 = "mmHg", paco2 = "mmHg", ph = "pH", lactate = "mmol/L", hb = "g/dL")
  allowed_states <- c("eligible_pending_candidate_build", "unavailable", "ambiguous",
                      "numeric_parse_hold", "source_scale_hold", "nonpositive_invalid",
                      "outside_broad_technical_range")
  for (name in ci1_metric_names()) {
    value_field <- paste0(name, "_candidate_value")
    unit_field <- paste0(name, "_canonical_unit")
    state_field <- paste0(name, "_value_state")
    x[[value_field]] <- ci1_number(x[[value_field]], code = "CI1_METRIC_NUMBER_INVALID")
    state <- x[[state_field]]
    if (any(is.na(state) | !state %in% allowed_states)) ci1_stop("CI1_METRIC_STATE_INVALID")
    observed <- state == "eligible_pending_candidate_build"
    if (any(observed & is.na(x[[value_field]])) || any(!observed & !is.na(x[[value_field]])))
      ci1_stop("CI1_METRIC_STATE_VALUE_MISMATCH")
    if (any(observed & (is.na(x[[unit_field]]) | x[[unit_field]] != units[[name]])))
      ci1_stop("CI1_METRIC_UNIT_INVALID")
    if (any(observed & (is.na(x[[paste0(name, "_linkage_class")]]) |
                        !nzchar(x[[paste0(name, "_linkage_class")]]))))
      ci1_stop("CI1_METRIC_LINKAGE_ABSENT")
    if (name == "ph" && any(observed & (x[[value_field]] < 6.3 | x[[value_field]] > 8.0)))
      ci1_stop("CI1_PH_RANGE_INVALID")
  }
  for (name in ci1_treatment_names()) {
    state <- x[[paste0(name, "_state")]]
    scope <- x[[paste0(name, "_scope")]]
    if (any(is.na(scope) | !scope %in% c("permitted_by_gate_r08", "excluded_by_gate_r08")))
      ci1_stop("CI1_TREATMENT_SCOPE_INVALID")
    if (any(scope == "permitted_by_gate_r08" & (is.na(state) | !state %in% c("yes", "unknown"))) ||
        any(scope == "excluded_by_gate_r08" & !is.na(state)))
      ci1_stop("CI1_TREATMENT_STATE_INVALID")
  }
  x
}

ci1_make_analytic_pairs <- function(matched_pairs, stays, clinical_pairs, cohort, window_minutes) {
  if (any(matched_pairs$cohort != cohort) || any(matched_pairs$pair_window_minutes != window_minutes))
    ci1_stop("CI1_PAIR_SCOPE_INVALID")
  pair_index <- match(ci1_key(matched_pairs, ci1_pair_keys()), ci1_key(clinical_pairs, ci1_pair_keys()))
  if (anyNA(pair_index)) ci1_stop("CI1_PAIR_ALIGNMENT_FAILED")
  side <- clinical_pairs[pair_index, , drop = FALSE]
  stay_index <- match(ci1_key(matched_pairs, c("patient_key_internal", "stay_key_internal")),
                      ci1_key(stays, c("patient_key_internal", "stay_key_internal")))
  if (anyNA(stay_index)) ci1_stop("CI1_STAY_ALIGNMENT_FAILED")
  stay <- stays[stay_index, setdiff(names(stays), c("patient_key_internal", "stay_key_internal")), drop = FALSE]
  out <- cbind(
    matched_pairs[c("cohort", "patient_key_internal", "stay_key_internal",
                    "spo2_stable_event_ordinal", "sao2_stable_event_ordinal",
                    "spo2_saturation_percent", "sao2_saturation_percent", "pair_window_minutes")],
    pair_key_internal = side$pair_key_internal,
    stay,
    side[setdiff(names(side), c(ci1_pair_keys(), "pair_key_internal"))]
  )
  out$bias_spo2_minus_sao2 <- out$spo2_saturation_percent - out$sao2_saturation_percent
  out$paired_mean_saturation_c90 <- (out$spo2_saturation_percent + out$sao2_saturation_percent) / 2 - 90
  if (nrow(out) != nrow(matched_pairs) || anyDuplicated(out$pair_key_internal))
    ci1_stop("CI1_ANALYTIC_CARDINALITY_INVALID")
  out
}
