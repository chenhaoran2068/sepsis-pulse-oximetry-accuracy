r211f_expected_cohorts <- function() c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")

r211f_protected_or_outcome_fields <- function() c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid", "bias_spo2_minus_sao2")

r211f_age_representation <- function(cohort) {
  out <- c(MIMIC = "numeric", Amsterdam = "interval", eICU = "numeric", SICDB = "interval", Lianyungang = "numeric")
  if (!cohort %in% names(out)) stop("R211F_COHORT_UNKNOWN", call. = FALSE)
  unname(out[[cohort]])
}

r211f_age_features <- function(cohort) {
  switch(cohort,
    MIMIC = data.frame(feature = c("age_value_per_10y", "age_topcoded_indicator"), feature_type = "continuous", feature_role = "candidate", stringsAsFactors = FALSE),
    eICU = data.frame(feature = c("age_value_per_10y", "age_topcoded_indicator"), feature_type = "continuous", feature_role = "candidate", stringsAsFactors = FALSE),
    Amsterdam = data.frame(feature = "age_interval", feature_type = "categorical", feature_role = "candidate", stringsAsFactors = FALSE),
    SICDB = data.frame(feature = "age_interval", feature_type = "categorical", feature_role = "candidate", stringsAsFactors = FALSE),
    Lianyungang = data.frame(feature = "age_value_per_10y", feature_type = "continuous", feature_role = "candidate", stringsAsFactors = FALSE),
    stop("R211F_COHORT_UNKNOWN", call. = FALSE)
  )
}

r211f_permitted_clinical_features <- function(cohort) {
  switch(cohort,
    MIMIC = c("concurrent_vasoactive_medication_use_state", "current_invasive_mechanical_ventilation_state"),
    Amsterdam = c("concurrent_vasoactive_medication_use_state", "current_invasive_mechanical_ventilation_state"),
    eICU = c("concurrent_vasoactive_medication_use_state"),
    SICDB = c("concurrent_vasoactive_medication_use_state"),
    Lianyungang = character(),
    stop("R211F_COHORT_UNKNOWN", call. = FALSE)
  )
}

r211f_feature_spec <- function(cohort) {
  age_row <- r211f_age_features(cohort)
  base <- rbind(
    data.frame(feature = "paired_mean_saturation_c90", feature_type = "continuous", feature_role = "structural", stringsAsFactors = FALSE),
    age_row,
    data.frame(
      feature = c("sex_class", "cci_value", "sofa_total", "pao2_candidate_value", "paco2_candidate_value", "ph_candidate_value", "lactate_candidate_value", "hb_candidate_value"),
      feature_type = c("categorical", rep("continuous", 7L)),
      feature_role = "candidate",
      stringsAsFactors = FALSE
    )
  )
  clinical <- r211f_permitted_clinical_features(cohort)
  if (length(clinical) > 0L) base <- rbind(base, data.frame(feature = clinical, feature_type = "categorical", feature_role = "candidate", stringsAsFactors = FALSE))
  rownames(base) <- NULL
  base
}

r211f_require_columns <- function(x, required, code) {
  if (!all(required %in% names(x))) stop(code, call. = FALSE)
  invisible(TRUE)
}

r211f_normalize_categorical <- function(value) {
  normalized <- tolower(trimws(as.character(value)))
  normalized[is.na(value) | normalized %in% c("", "unknown", "unknown_or_missing", "na", "n/a")] <- NA_character_
  normalized
}

r211f_as_continuous <- function(value) {
  raw <- trimws(as.character(value))
  absent <- is.na(value) | raw == ""
  parsed <- suppressWarnings(as.numeric(value))
  if (any(!absent & is.na(parsed))) stop("R211F_CONTINUOUS_NONNUMERIC", call. = FALSE)
  if (any(!is.na(parsed) & !is.finite(parsed))) stop("R211F_CONTINUOUS_NONFINITE", call. = FALSE)
  parsed
}

r211f_validate_candidate <- function(candidate, cohort) {
  spec <- r211f_feature_spec(cohort)
  required <- unique(c(r211f_protected_or_outcome_fields(), "spo2_percent", "sao2_percent", spec$feature))
  r211f_require_columns(candidate, required, "R211F_CANDIDATE_INTERFACE_MISSING")
  if (nrow(candidate) == 0L || anyNA(candidate$cohort) || !all(as.character(candidate$cohort) == cohort)) stop("R211F_COHORT_IDENTITY_INVALID", call. = FALSE)
  if (anyNA(candidate$analysis_patient_key) || anyNA(candidate$analysis_stay_key) || anyNA(candidate$pair_uid) || anyDuplicated(candidate$pair_uid)) stop("R211F_PROTECTED_KEY_INTEGRITY_FAILURE", call. = FALSE)
  if (any(!is.finite(candidate$spo2_percent)) || any(!is.finite(candidate$sao2_percent)) || any(!is.finite(candidate$bias_spo2_minus_sao2))) stop("R211F_STRUCTURAL_OXYGEN_OR_BIAS_MISSING", call. = FALSE)
  if (!isTRUE(all.equal(candidate$bias_spo2_minus_sao2, candidate$spo2_percent - candidate$sao2_percent, check.attributes = FALSE))) stop("R211F_BIAS_DIRECTION_INVALID", call. = FALSE)
  if (any(!is.finite(candidate$paired_mean_saturation_c90)) || !isTRUE(all.equal(candidate$paired_mean_saturation_c90, ((candidate$spo2_percent + candidate$sao2_percent) / 2) - 90, check.attributes = FALSE))) stop("R211F_PAIRED_MEAN_IDENTITY_INVALID", call. = FALSE)
  if (cohort %in% c("MIMIC", "eICU") && any(!candidate$age_topcoded_indicator %in% c(0, 1), na.rm = TRUE)) stop("R211F_AGE_TOPCODE_INDICATOR_INVALID", call. = FALSE)
  invisible(TRUE)
}

r211f_observed <- function(value, feature_type) {
  if (feature_type == "continuous") return(!is.na(r211f_as_continuous(value)))
  !is.na(r211f_normalize_categorical(value))
}

r211f_mode <- function(value) {
  observed <- value[!is.na(value)]
  if (length(observed) == 0L) stop("R211F_CATEGORICAL_MODE_UNAVAILABLE", call. = FALSE)
  tab <- table(observed)
  sort(names(tab)[tab == max(tab)])[[1L]]
}

r211f_validate_support_parameters <- function(minimum_support, minimum_observed_n) {
  if (!isTRUE(all.equal(minimum_support, 0.10)) || length(minimum_support) != 1L) stop("R211F_SUPPORT_PROPORTION_NOT_APPROVED", call. = FALSE)
  if (!identical(as.integer(minimum_observed_n), 500L) || length(minimum_observed_n) != 1L) stop("R211F_SUPPORT_COUNT_NOT_APPROVED", call. = FALSE)
  invisible(TRUE)
}

r211f_fit_training_recipe <- function(training_data, cohort, minimum_support = 0.10, minimum_observed_n = 500L) {
  r211f_validate_support_parameters(minimum_support, minimum_observed_n)
  r211f_validate_candidate(training_data, cohort)
  spec <- r211f_feature_spec(cohort)
  rows <- lapply(seq_len(nrow(spec)), function(i) {
    feature <- spec$feature[[i]]
    type <- spec$feature_type[[i]]
    role <- spec$feature_role[[i]]
    if (role == "structural") {
      return(data.frame(feature = feature, feature_type = type, feature_role = role, observed_n = nrow(training_data), observed_support = 1, retained = TRUE, fill_value = NA_character_, allowed_levels = NA_character_, stringsAsFactors = FALSE))
    }
    raw_value <- training_data[[feature]]
    value <- if (type == "continuous") r211f_as_continuous(raw_value) else r211f_normalize_categorical(raw_value)
    observed <- !is.na(value)
    observed_n <- sum(observed)
    observed_support <- observed_n / nrow(training_data)
    retained <- observed_support >= minimum_support && observed_n >= minimum_observed_n
    fill_value <- NA_character_
    allowed_levels <- NA_character_
    if (retained && type == "continuous") fill_value <- as.character(stats::median(value[observed]))
    if (retained && type == "categorical") {
      fill_value <- r211f_mode(value)
      allowed_levels <- paste(sort(unique(value[observed])), collapse = "|")
    }
    data.frame(feature = feature, feature_type = type, feature_role = role, observed_n = observed_n, observed_support = observed_support, retained = retained, fill_value = fill_value, allowed_levels = allowed_levels, stringsAsFactors = FALSE)
  })
  recipe <- do.call(rbind, rows)
  attr(recipe, "cohort") <- cohort
  attr(recipe, "minimum_support") <- minimum_support
  attr(recipe, "minimum_observed_n") <- as.integer(minimum_observed_n)
  recipe
}

r211f_validate_recipe <- function(recipe, cohort) {
  r211f_require_columns(recipe, c("feature", "feature_type", "feature_role", "observed_n", "observed_support", "retained", "fill_value", "allowed_levels"), "R211F_RECIPE_INTERFACE_MISSING")
  spec <- r211f_feature_spec(cohort)
  if (nrow(recipe) != nrow(spec) || !setequal(recipe$feature, spec$feature) || any(recipe$feature %in% r211f_protected_or_outcome_fields())) stop("R211F_RECIPE_FEATURE_SCOPE_INVALID", call. = FALSE)
  if (!identical(attr(recipe, "cohort"), cohort) || !isTRUE(all.equal(attr(recipe, "minimum_support"), 0.10)) || !identical(as.integer(attr(recipe, "minimum_observed_n")), 500L)) stop("R211F_RECIPE_PARAMETER_INVALID", call. = FALSE)
  structural <- recipe$feature_role == "structural"
  if (!all(recipe$retained[structural]) || any(!is.na(recipe$fill_value[structural]))) stop("R211F_STRUCTURAL_RECIPE_INVALID", call. = FALSE)
  if (any(recipe$retained & recipe$feature_role == "candidate" & (recipe$observed_support < 0.10 | recipe$observed_n < 500L))) stop("R211F_SUPPORT_RULE_BYPASSED", call. = FALSE)
  invisible(TRUE)
}

r211f_apply_fixed_recipe <- function(data, recipe, cohort) {
  r211f_validate_candidate(data, cohort)
  r211f_validate_recipe(recipe, cohort)
  output <- data[, c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid", "bias_spo2_minus_sao2"), drop = FALSE]
  for (i in seq_len(nrow(recipe))) {
    feature <- recipe$feature[[i]]
    type <- recipe$feature_type[[i]]
    role <- recipe$feature_role[[i]]
    if (role == "structural") {
      output[[paste0(feature, "__model")]] <- data[[feature]]
      next
    }
    if (!isTRUE(recipe$retained[[i]])) next
    value <- if (type == "continuous") r211f_as_continuous(data[[feature]]) else r211f_normalize_categorical(data[[feature]])
    missing <- is.na(value)
    output[[paste0(feature, "__was_missing")]] <- missing
    if (type == "continuous") {
      model_value <- value
      model_value[missing] <- as.numeric(recipe$fill_value[[i]])
    } else {
      allowed <- strsplit(recipe$allowed_levels[[i]], "|", fixed = TRUE)[[1L]]
      unseen <- !missing & !value %in% allowed
      model_value <- value
      model_value[missing] <- recipe$fill_value[[i]]
      model_value[unseen] <- "__unseen__"
      output[[paste0(feature, "__was_unseen")]] <- unseen
    }
    output[[paste0(feature, "__model")]] <- model_value
  }
  output
}
