rfp_stop_if <- function(condition, code) {
  if (isTRUE(condition)) stop(code, call. = FALSE)
  invisible(TRUE)
}

rfp_round_c90_to_whole_percentage_point <- function(value) {
  value <- as.numeric(value)
  rfp_stop_if(anyNA(value) || any(!is.finite(value)), "RFP_PAIRED_MEAN_INVALID")
  rounded <- floor((value + 90) + 0.5) - 90
  rfp_stop_if(any(abs(rounded - round(rounded)) > 1e-12), "RFP_ROUNDED_VALUE_NONINTEGER")
  rounded
}

rfp_transform_model_frame <- function(frame) {
  field <- "paired_mean_saturation_c90__model"
  rfp_stop_if(!field %in% names(frame), "RFP_MODEL_PAIRED_MEAN_ABSENT")
  outcome_before <- frame$bias_spo2_minus_sao2
  other_fields <- setdiff(names(frame), field)
  other_before <- frame[, other_fields, drop = FALSE]
  frame[[field]] <- rfp_round_c90_to_whole_percentage_point(frame[[field]])
  rfp_stop_if(!identical(outcome_before, frame$bias_spo2_minus_sao2), "RFP_OUTCOME_CHANGED")
  rfp_stop_if(!identical(other_before, frame[, other_fields, drop = FALSE]), "RFP_NONSTRUCTURAL_FIELD_CHANGED")
  frame
}

rfp_transform_source_features <- function(frame) {
  field <- "paired_mean_saturation_c90"
  rfp_stop_if(!field %in% names(frame), "RFP_SOURCE_PAIRED_MEAN_ABSENT")
  frame[[field]] <- rfp_round_c90_to_whole_percentage_point(frame[[field]])
  frame
}

rfp_precision_diagnostic <- function(frame, cohort, partition, fold = NA_integer_) {
  field <- "paired_mean_saturation_c90__model"
  rfp_stop_if(!field %in% names(frame), "RFP_DIAGNOSTIC_FIELD_ABSENT")
  value <- as.numeric(frame[[field]])
  data.frame(
    cohort = cohort,
    partition = partition,
    inner_validation_fold = as.integer(fold),
    row_n = length(value),
    unique_level_n = length(unique(value)),
    noninteger_n = sum(abs(value - round(value)) > 1e-12),
    minimum = min(value),
    maximum = max(value),
    representation = "nearest_whole_percentage_point_centered_at_90",
    stringsAsFactors = FALSE
  )
}
