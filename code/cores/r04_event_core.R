r04_stop <- function(code) {
  stop(structure(list(message = code, call = NULL), class = c("r04_error", "error", "condition")))
}

r04_text <- function(x) {
  out <- trimws(as.character(x))
  out[is.na(out)] <- ""
  out
}

r04_num <- function(x) suppressWarnings(as.numeric(r04_text(x)))

r04_sha256_text <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) r04_stop("R04_DIGEST_DEPENDENCY_ABSENT")
  vapply(enc2utf8(as.character(x)), digest::digest, character(1L), algo = "sha256", serialize = FALSE)
}

r04_internal_key <- function(cohort, kind, raw_key) {
  paste0(tolower(cohort), "_", kind, "_", r04_sha256_text(paste(cohort, kind, r04_text(raw_key), sep = "|")))
}

r04_named_lookup_or_null <- function(x, name) {
  out <- unname(x[name])
  if (!length(out) || is.na(out[[1L]]) || !nzchar(as.character(out[[1L]]))) return(NULL)
  as.character(out[[1L]])
}

r04_range_status <- function(x) {
  value <- r04_num(x)
  out <- rep("not_finite", length(value))
  out[is.finite(value) & value < 70] <- "below_70"
  out[is.finite(value) & value >= 70 & value <= 100] <- "eligible_70_100"
  out[is.finite(value) & value > 100] <- "above_100"
  out
}

r04_unit_percent <- function(x) {
  tolower(r04_text(x)) %in% c("%", "percent", "pct")
}

r04_mimic_asset_exception_counts <- function(asset_exception, time_ok) {
  if (length(asset_exception) != length(time_ok)) r04_stop("R04_MIMIC_ASSET_COUNT_LENGTH")
  c(
    all_p1_source_rows = sum(asset_exception, na.rm = TRUE),
    inside_icu_interval = sum(asset_exception & time_ok, na.rm = TRUE)
  )
}

r04_mimic_sao2_numeric_scope_counts <- function(numeric_bad, time_ok) {
  if (length(numeric_bad) != length(time_ok)) r04_stop("R04_MIMIC_SAO2_COUNT_LENGTH")
  if (anyNA(numeric_bad) || anyNA(time_ok)) r04_stop("R04_MIMIC_SAO2_COUNT_NA")
  c(
    p1_hospital_rows = length(numeric_bad),
    p1_hospital_missing_or_nonnumeric = sum(numeric_bad),
    p1_icu_missing_or_nonnumeric = sum(numeric_bad & time_ok),
    p1_outside_or_invalid_icu_missing_or_nonnumeric = sum(numeric_bad & !time_ok)
  )
}

r04_qa_rows <- function(cohort, modality, entries) {
  if (is.null(names(entries)) || any(!nzchar(names(entries))) || any(lengths(entries) != 1L)) r04_stop("R04_QA_NAMED_SCALAR_REQUIRED")
  data.frame(
    cohort = cohort,
    modality = modality,
    metric = names(entries),
    n = as.double(unlist(entries, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

r04_mimic_subject_hadm_text_keys <- function(x) {
  required <- c("subject_id", "hadm_id")
  if (!is.data.frame(x) || length(setdiff(required, names(x)))) r04_stop("R04_MIMIC_SUBJECT_HADM_SCHEMA")
  x$subject_id <- r04_text(x$subject_id)
  x$hadm_id <- r04_text(x$hadm_id)
  x
}

r04_parse_utc_min <- function(x) {
  if (inherits(x, "POSIXt")) return(as.numeric(x) / 60)
  as.numeric(as.POSIXct(r04_text(x), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) / 60
}

r04_parse_lyy_utc_min <- function(x) {
  if (!requireNamespace("lubridate", quietly = TRUE)) r04_stop("R04_LUBRIDATE_DEPENDENCY_ABSENT")
  as.numeric(lubridate::parse_date_time(r04_text(x), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC", quiet = TRUE)) / 60
}

r04_parse_lyy_value <- function(x) {
  if (!requireNamespace("stringr", quietly = TRUE)) r04_stop("R04_STRINGR_DEPENDENCY_ABSENT")
  suppressWarnings(as.numeric(stringr::str_extract(as.character(x), "[0-9.]+")))
}

r04_mimic_label_status <- function(x) {
  y <- toupper(r04_text(x))
  ven <- grepl("VEN", y, fixed = TRUE)
  mix <- grepl("MIX", y, fixed = TRUE)
  art <- grepl("ART", y, fixed = TRUE)
  out <- rep("unclassified_nonempty", length(y))
  out[!nzchar(y)] <- "unlabelled_operational_candidate"
  out[art & !ven & !mix] <- "arterial_explicit"
  out[ven & !mix] <- "venous_explicit"
  out[mix & !ven] <- "mixed_explicit"
  out[ven & mix] <- "ambiguous_or_conflict"
  out
}

r04_mimic_specimen_lookup <- function(labels) {
  if (!requireNamespace("data.table", quietly = TRUE)) r04_stop("R04_DATATABLE_DEPENDENCY_ABSENT")
  required <- c("subject_id", "hadm_id", "specimen_id", "value")
  if (!is.data.frame(labels) || length(setdiff(required, names(labels)))) r04_stop("R04_MIMIC_LABEL_SCHEMA")
  x <- data.table::as.data.table(labels)
  x[, specimen_link_key := paste(r04_text(subject_id), r04_text(hadm_id), r04_text(specimen_id), sep = "|")]
  x[!nzchar(r04_text(subject_id)) | !nzchar(r04_text(hadm_id)) | !nzchar(r04_text(specimen_id)), specimen_link_key := NA_character_]
  x[, label_status := r04_mimic_label_status(value)]
  x <- x[!is.na(specimen_link_key)]
  if (!nrow(x)) return(data.frame(specimen_link_key = character(), specimen_status = character(), stringsAsFactors = FALSE))
  x[, specimen_status := if (data.table::uniqueN(label_status) == 1L) label_status[[1L]] else "ambiguous_or_conflict", by = specimen_link_key]
  unique(x[, .(specimen_link_key, specimen_status)])
}

r04_amsterdam_is_venous_comment <- function(x) {
  y <- tolower(r04_text(x))
  grepl("venous", y, fixed = TRUE) | grepl("veneus", y, fixed = TRUE)
}

r04_amsterdam_normalize_sao2 <- function(itemid, value) {
  item <- suppressWarnings(as.integer(r04_text(itemid)))
  out <- r04_num(value)
  convert <- item == 12311L & is.finite(out) & out > 0 & out <= 1
  out[convert] <- out[convert] * 100
  list(value = out, fraction_normalized = convert)
}

r04_lyy_spo2_label_keep <- function(x) {
  item <- r04_text(x)
  grepl("血氧饱和度", item, fixed = TRUE) | grepl("SpO2", item, fixed = TRUE)
}

r04_lyy_sao2_label_keep <- function(x) {
  item <- r04_text(x)
  item == "氧饱和度(SaO2)-动脉血"
}

r04_build_p1_bridge <- function(cohort, p1) {
  required <- c("stay_key", "patient_key")
  if (!is.data.frame(p1) || length(setdiff(required, names(p1)))) r04_stop("R04_P1_SCHEMA")
  out <- data.frame(
    stay_raw = r04_text(p1$stay_key),
    patient_raw = r04_text(p1$patient_key),
    stringsAsFactors = FALSE
  )
  if (!nrow(out) || any(!nzchar(out$stay_raw)) || any(!nzchar(out$patient_raw)) || anyDuplicated(out$stay_raw)) r04_stop("R04_P1_KEY_INTERFACE")
  out$stay_key_internal <- r04_internal_key(cohort, "stay", out$stay_raw)
  out$patient_key_internal <- r04_internal_key(cohort, "patient", out$patient_raw)
  out
}

r04_finalize_events <- function(x, cohort, modality, artifact_threshold = NA_real_, artifact_assessable = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) r04_stop("R04_DATATABLE_DEPENDENCY_ABSENT")
  required <- c("stay_raw", "stay_key_internal", "patient_key_internal", "event_time_min", "saturation_percent", "source_family", "source_priority")
  if (!is.data.frame(x) || length(setdiff(required, names(x)))) r04_stop("R04_EVENT_STAGE_SCHEMA")
  d <- data.table::copy(data.table::as.data.table(x))
  d[, stay_raw := r04_text(stay_raw)]
  d[, event_time_min := r04_num(event_time_min)]
  d[, saturation_percent := r04_num(saturation_percent)]
  d[, source_priority := suppressWarnings(as.integer(source_priority))]
  d <- d[nzchar(stay_raw) & is.finite(event_time_min) & is.finite(saturation_percent) & !is.na(source_priority)]
  d[, range_status := r04_range_status(saturation_percent)]
  d[, analysis_eligible_range := range_status == "eligible_70_100"]
  input_valid <- nrow(d)
  data.table::setorder(d, stay_raw, event_time_min, saturation_percent, source_priority, source_family)
  d[, exact_duplicate_source_signature := duplicated(d, by = c("stay_raw", "event_time_min", "saturation_percent", "source_family"))]
  d <- d[exact_duplicate_source_signature == FALSE]
  d[, cross_source_exact_duplicate_drop := FALSE]
  if (modality == "SaO2" && nrow(d)) {
    d[, min_source_priority := min(source_priority), by = .(stay_raw, event_time_min, saturation_percent)]
    d[, cross_source_exact_duplicate_drop := source_priority > min_source_priority]
    d <- d[cross_source_exact_duplicate_drop == FALSE]
    d[, min_source_priority := NULL]
  }
  data.table::setorder(d, stay_raw, event_time_min, source_priority, saturation_percent, source_family)
  d[, stable_event_ordinal := seq_len(.N), by = stay_raw]
  d[, artifact_qc_status := "not_applicable_non_spo2"]
  d[, possible_transient_artifact := FALSE]
  if (identical(modality, "SpO2")) {
    if (isTRUE(artifact_assessable)) {
      if (!is.finite(artifact_threshold) || artifact_threshold <= 0) r04_stop("R04_ARTIFACT_THRESHOLD_INTERFACE")
      d[, `:=`(artifact_qc_status = "assessed_not_flagged", possible_transient_artifact = FALSE)]
      eligible <- d$analysis_eligible_range
      if (any(eligible)) {
        z <- d[eligible]
        data.table::setorder(z, stay_raw, event_time_min, stable_event_ordinal)
        z[, `:=`(
          prev_time = data.table::shift(event_time_min),
          prev_value = data.table::shift(saturation_percent),
          next_time = data.table::shift(event_time_min, type = "lead"),
          next_value = data.table::shift(saturation_percent, type = "lead")
        ), by = stay_raw]
        z[, residual := abs(saturation_percent - (prev_value + (next_value - prev_value) * ((event_time_min - prev_time) / (next_time - prev_time))))]
        z[, artifact_flag := is.finite(prev_time) & is.finite(next_time) & (event_time_min - prev_time) > 0 & (next_time - event_time_min) > 0 &
          (event_time_min - prev_time) <= 15 & (next_time - event_time_min) <= 15 &
          (saturation_percent > pmax(prev_value, next_value) | saturation_percent < pmin(prev_value, next_value)) &
          is.finite(residual) & residual >= artifact_threshold]
        d[eligible, possible_transient_artifact := z$artifact_flag]
        d[eligible & possible_transient_artifact, artifact_qc_status := "possible_transient_artifact"]
      }
    } else {
      d[, artifact_qc_status := "not_assessable_temporal_resolution"]
    }
  }
  out <- d[, .(stay_key_internal, patient_key_internal, event_time_min, saturation_percent, source_family, source_priority, range_status, analysis_eligible_range, artifact_qc_status, possible_transient_artifact, stable_event_ordinal)]
  if (anyDuplicated(out, by = c("stay_key_internal", "event_time_min", "saturation_percent", "source_family", "source_priority", "stable_event_ordinal"))) r04_stop("R04_FINAL_EVENT_DUPLICATE")
  list(
    events = out,
    aggregate = data.frame(
      cohort = cohort,
      modality = modality,
      input_valid_event_rows = as.double(input_valid),
      exact_duplicate_rows_removed = as.double(input_valid - nrow(d)),
      final_candidate_events = as.double(nrow(out)),
      below_70_events = as.double(sum(out$range_status == "below_70")),
      eligible_70_100_events = as.double(sum(out$range_status == "eligible_70_100")),
      above_100_events = as.double(sum(out$range_status == "above_100")),
      possible_transient_artifact_events = as.double(sum(out$possible_transient_artifact)),
      artifact_qc_mode = if (identical(modality, "SpO2")) if (isTRUE(artifact_assessable)) "P99_15_MIN_THREE_POINT" else "NOT_ASSESSABLE_TEMPORAL_RESOLUTION" else "NOT_APPLICABLE_SAO2",
      stringsAsFactors = FALSE
    )
  )
}

r04_assert_controlled_event_output <- function(events) {
  required <- c("stay_key_internal", "patient_key_internal", "event_time_min", "saturation_percent", "source_family", "source_priority", "range_status", "analysis_eligible_range", "artifact_qc_status", "possible_transient_artifact", "stable_event_ordinal")
  if (!is.data.frame(events) || length(setdiff(required, names(events)))) r04_stop("R04_CONTROLLED_EVENT_SCHEMA")
  forbidden <- c("stay_raw", "patient_raw", "subject_id", "hadm_id", "stay_id", "admissionid", "patientunitstayid", "CaseID", "patient_sn", "event_time_absolute", "comment", "specimen_label")
  if (length(intersect(forbidden, names(events)))) r04_stop("R04_CONTROLLED_EVENT_RAW_FIELD")
  if (any(!nzchar(r04_text(events$stay_key_internal))) || any(!nzchar(r04_text(events$patient_key_internal)))) r04_stop("R04_CONTROLLED_EVENT_KEY_BLANK")
  if (any(!is.finite(r04_num(events$event_time_min))) || any(!is.finite(r04_num(events$saturation_percent)))) r04_stop("R04_CONTROLLED_EVENT_VALUE_INVALID")
  invisible(TRUE)
}
