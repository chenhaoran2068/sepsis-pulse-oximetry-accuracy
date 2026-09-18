si_stop <- function(code) stop(code, call. = FALSE)

si_required_stay_columns <- c("stay_key_internal", "patient_key_internal")
si_required_event_columns <- c(
  "stay_key_internal", "patient_key_internal", "event_time_min", "saturation_percent",
  "source_family", "source_priority", "range_status", "analysis_eligible_range",
  "artifact_qc_status", "possible_transient_artifact", "stable_event_ordinal"
)

si_read_tsv <- function(path, expected_columns, kind) {
  if (!file.exists(path) || dir.exists(path)) si_stop(paste0("SI_", kind, "_FILE_MISSING"))
  out <- tryCatch(utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
    comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    fileEncoding = "UTF-8"), error = function(e) si_stop(paste0("SI_", kind, "_READ_ERROR")))
  if (!identical(names(out), expected_columns)) si_stop(paste0("SI_", kind, "_SCHEMA_ERROR"))
  if (!nrow(out)) si_stop(paste0("SI_", kind, "_EMPTY"))
  out
}

si_assert_text <- function(x, code) {
  if (anyNA(x) || any(!nzchar(trimws(x))) || any(trimws(x) != x)) si_stop(code)
  invisible(TRUE)
}

si_parse_numeric <- function(x, code, integer = FALSE) {
  if (anyNA(x) || any(!nzchar(x))) si_stop(code)
  parsed <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(parsed))) si_stop(code)
  if (integer && any(parsed != floor(parsed))) si_stop(code)
  parsed
}

si_parse_logical <- function(x, code) {
  if (anyNA(x) || any(!x %in% c("TRUE", "FALSE"))) si_stop(code)
  x == "TRUE"
}

si_validate_stays <- function(stays) {
  if (!is.data.frame(stays) || !identical(names(stays), si_required_stay_columns)) si_stop("SI_STAYS_SCHEMA_ERROR")
  if (!nrow(stays)) si_stop("SI_STAYS_EMPTY")
  si_assert_text(stays$stay_key_internal, "SI_STAY_KEY_INVALID")
  si_assert_text(stays$patient_key_internal, "SI_PATIENT_KEY_INVALID")
  if (anyDuplicated(stays$stay_key_internal)) si_stop("SI_STAY_KEY_DUPLICATE")
  stays
}

si_validate_events <- function(events, stays, cohort, modality) {
  if (!is.data.frame(events) || !identical(names(events), si_required_event_columns)) si_stop("SI_EVENT_SCHEMA_ERROR")
  if (!nrow(events)) si_stop(paste0("SI_", modality, "_EMPTY"))
  si_assert_text(events$stay_key_internal, "SI_EVENT_STAY_KEY_INVALID")
  si_assert_text(events$patient_key_internal, "SI_EVENT_PATIENT_KEY_INVALID")
  si_assert_text(events$source_family, "SI_SOURCE_FAMILY_INVALID")
  si_assert_text(events$range_status, "SI_RANGE_STATUS_INVALID")
  si_assert_text(events$artifact_qc_status, "SI_ARTIFACT_STATUS_INVALID")
  index <- match(events$stay_key_internal, stays$stay_key_internal)
  if (anyNA(index)) si_stop("SI_EVENT_STAY_NOT_IN_ROSTER")
  if (any(events$patient_key_internal != stays$patient_key_internal[index])) si_stop("SI_EVENT_PATIENT_MISMATCH")
  events$event_time_min <- si_parse_numeric(events$event_time_min, "SI_EVENT_TIME_INVALID")
  if (any(events$event_time_min < 0)) si_stop("SI_EVENT_TIME_INVALID")
  events$saturation_percent <- si_parse_numeric(events$saturation_percent, "SI_SATURATION_INVALID")
  events$source_priority <- si_parse_numeric(events$source_priority, "SI_SOURCE_PRIORITY_INVALID", integer = TRUE)
  if (any(events$source_priority < 1)) si_stop("SI_SOURCE_PRIORITY_INVALID")
  events$stable_event_ordinal <- si_parse_numeric(events$stable_event_ordinal, "SI_EVENT_ORDINAL_INVALID", integer = TRUE)
  events$analysis_eligible_range <- si_parse_logical(events$analysis_eligible_range, "SI_ELIGIBLE_FLAG_INVALID")
  events$possible_transient_artifact <- si_parse_logical(events$possible_transient_artifact, "SI_ARTIFACT_FLAG_INVALID")
  if (any(events$range_status != r04_range_status(events$saturation_percent))) si_stop("SI_RANGE_VALUE_MISMATCH")
  if (modality == "SaO2" && (any(events$possible_transient_artifact) ||
      any(events$artifact_qc_status != "not_applicable_non_spo2"))) si_stop("SI_SAO2_ARTIFACT_STATE_INVALID")
  if (modality == "SpO2" && any(events$possible_transient_artifact !=
      (events$artifact_qc_status == "possible_transient_artifact"))) si_stop("SI_SPO2_ARTIFACT_STATE_INVALID")
  r05_validate_events(events, cohort, modality)
  events
}

si_load_input <- function(input_dir, cohort) {
  if (length(cohort) != 1L || !cohort %in% c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")) si_stop("SI_COHORT_INVALID")
  if (!dir.exists(input_dir)) si_stop("SI_INPUT_DIRECTORY_MISSING")
  stays <- si_validate_stays(si_read_tsv(file.path(input_dir, "stays.tsv"), si_required_stay_columns, "STAYS"))
  spo2 <- si_validate_events(si_read_tsv(file.path(input_dir, "spo2_events.tsv"), si_required_event_columns, "SPO2"), stays, cohort, "SpO2")
  sao2 <- si_validate_events(si_read_tsv(file.path(input_dir, "sao2_events.tsv"), si_required_event_columns, "SAO2"), stays, cohort, "SaO2")
  r05_validate_modalities(spo2, sao2)
  list(stays = stays, spo2 = spo2, sao2 = sao2)
}
