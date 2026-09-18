rfag_stop_if <- function(condition, code) {
  if (isTRUE(condition)) stop(code, call. = FALSE)
  invisible(TRUE)
}

rfag_hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

rfag_protected_fields <- function() c(
  "cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid",
  "bias_spo2_minus_sao2", "ml_input_partition", "inner_validation_fold"
)

rfag_forbidden_oxygen_fields <- function() c(
  "spo2_percent", "sao2_percent", "spo2_percent__model", "sao2_percent__model"
)

rfag_predictors <- function(x) {
  predictors <- grep("__(model|was_missing|was_unseen)$", names(x), value = TRUE)
  predictors <- setdiff(predictors, rfag_protected_fields())
  rfag_stop_if(length(predictors) == 0L, "RFAG_PREDICTOR_SET_EMPTY")
  rfag_stop_if(!"paired_mean_saturation_c90__model" %in% predictors, "RFAG_PAIRED_MEAN_ABSENT")
  rfag_stop_if(any(rfag_forbidden_oxygen_fields() %in% predictors), "RFAG_STANDALONE_OXYGEN_PREDICTOR_PRESENT")
  predictors
}

rfag_validate_frame <- function(x, cohort, partition, fold) {
  required <- rfag_protected_fields()
  rfag_stop_if(!all(required %in% names(x)), "RFAG_INPUT_INTERFACE_MISSING")
  rfag_stop_if(any(rfag_forbidden_oxygen_fields() %in% names(x)), "RFAG_FORBIDDEN_OXYGEN_FIELD_PRESENT")
  rfag_stop_if(nrow(x) == 0L, "RFAG_INPUT_EMPTY")
  rfag_stop_if(anyNA(x$analysis_patient_key) || anyNA(x$analysis_stay_key) || anyNA(x$pair_uid), "RFAG_KEY_MISSING")
  rfag_stop_if(anyDuplicated(x$pair_uid) > 0L, "RFAG_PAIR_DUPLICATED")
  rfag_stop_if(!all(as.character(x$cohort) == cohort), "RFAG_COHORT_INVALID")
  rfag_stop_if(!all(as.character(x$ml_input_partition) == partition), "RFAG_PARTITION_INVALID")
  rfag_stop_if(!all(as.integer(x$inner_validation_fold) == as.integer(fold)), "RFAG_FOLD_INVALID")
  rfag_stop_if(any(!is.finite(x$bias_spo2_minus_sao2)), "RFAG_OUTCOME_INVALID")
  predictors <- rfag_predictors(x)
  for (column in predictors) {
    value <- x[[column]]
    if (is.numeric(value) || is.integer(value) || is.logical(value)) {
      rfag_stop_if(anyNA(value) || any(!is.finite(as.numeric(value))), paste0("RFAG_NUMERIC_INVALID:", column))
    } else if (is.character(value) || is.factor(value)) {
      normalized <- trimws(as.character(value))
      rfag_stop_if(anyNA(value) || any(normalized == ""), paste0("RFAG_CATEGORICAL_INVALID:", column))
    } else {
      stop(paste0("RFAG_TYPE_UNSUPPORTED:", column), call. = FALSE)
    }
  }
  invisible(predictors)
}

rfag_recipe_for_fold <- function(recipe_bundle, cohort, fold) {
  required <- c("recipe_context", "feature", "feature_type", "feature_role", "retained", "allowed_levels")
  rfag_stop_if(!all(required %in% names(recipe_bundle)), "RFAG_RECIPE_INTERFACE_MISSING")
  context <- paste0("inner_fold_", as.integer(fold))
  recipe <- recipe_bundle[as.character(recipe_bundle$recipe_context) == context, , drop = FALSE]
  rfag_stop_if(nrow(recipe) == 0L || anyDuplicated(recipe$feature) > 0L, "RFAG_RECIPE_CONTEXT_INVALID")
  rfag_stop_if(!any(recipe$feature == "paired_mean_saturation_c90" & recipe$feature_role == "structural" & recipe$retained), "RFAG_RECIPE_PAIRED_MEAN_INVALID")
  recipe
}

rfag_factor_spec <- function(recipe) {
  result <- list()
  rows <- recipe$feature_type == "categorical" & recipe$retained
  for (index in which(rows)) {
    feature <- as.character(recipe$feature[[index]])
    allowed <- strsplit(as.character(recipe$allowed_levels[[index]]), "|", fixed = TRUE)[[1L]]
    allowed <- unique(c(allowed, "__unseen__"))
    result[[paste0(feature, "__model")]] <- list(
      levels = allowed,
      ordered = identical(feature, "age_interval")
    )
  }
  result
}

rfag_prepare_frame <- function(x, predictors, factor_spec) {
  output <- x[, c("bias_spo2_minus_sao2", predictors), drop = FALSE]
  for (column in predictors) {
    value <- output[[column]]
    if (column %in% names(factor_spec)) {
      spec <- factor_spec[[column]]
      value <- trimws(as.character(value))
      rfag_stop_if(any(!value %in% spec$levels), paste0("RFAG_LEVEL_OUTSIDE_RECIPE:", column))
      output[[column]] <- factor(value, levels = spec$levels, ordered = isTRUE(spec$ordered))
    } else {
      if (is.logical(value)) value <- as.integer(value)
      output[[column]] <- as.numeric(value)
      rfag_stop_if(any(!is.finite(output[[column]])), paste0("RFAG_PREPARED_NUMERIC_INVALID:", column))
    }
  }
  output
}

rfag_prepare_fold <- function(training, validation, recipe, cohort, fold) {
  train_predictors <- rfag_validate_frame(training, cohort, "inner_training", fold)
  valid_predictors <- rfag_validate_frame(validation, cohort, "inner_validation", fold)
  rfag_stop_if(!identical(train_predictors, valid_predictors), "RFAG_PREDICTOR_ORDER_MISMATCH")
  rfag_stop_if(length(intersect(as.character(training$analysis_patient_key), as.character(validation$analysis_patient_key))) > 0L, "RFAG_PATIENT_LEAKAGE")
  rfag_stop_if(length(intersect(as.character(training$pair_uid), as.character(validation$pair_uid))) > 0L, "RFAG_PAIR_LEAKAGE")
  spec <- rfag_factor_spec(recipe)
  list(
    predictors = train_predictors,
    factor_spec = spec,
    training = rfag_prepare_frame(training, train_predictors, spec),
    validation = rfag_prepare_frame(validation, train_predictors, spec)
  )
}

rfag_resolve_mtry <- function(p, rule = "p_over_3") {
  divisors <- c(p_over_4 = 4, p_over_3 = 3, p_over_2 = 2, p = 1)
  rfag_stop_if(length(p) != 1L || !is.finite(p) || p < 1L || !rule %in% names(divisors), "RFAG_MTRY_INVALID")
  as.integer(max(1L, min(as.integer(p), round(as.integer(p) / divisors[[rule]]))))
}

rfag_fit_predict <- function(training, validation, mtry, min_node_size, seed, num_trees = 500L, num_threads = 1L) {
  rfag_stop_if(!requireNamespace("ranger", quietly = TRUE), "RFAG_RANGER_UNAVAILABLE")
  rfag_stop_if(!num_threads %in% c(1L, 4L), "RFAG_THREAD_COUNT_INVALID")
  started <- proc.time()[[3L]]
  fit <- ranger::ranger(
    dependent.variable.name = "bias_spo2_minus_sao2",
    data = training,
    num.trees = as.integer(num_trees),
    mtry = as.integer(mtry),
    min.node.size = as.integer(min_node_size),
    splitrule = "variance",
    replace = TRUE,
    respect.unordered.factors = "partition",
    importance = "none",
    write.forest = TRUE,
    num.threads = as.integer(num_threads),
    seed = as.integer(seed)
  )
  prediction <- as.numeric(predict(fit, data = validation, num.threads = as.integer(num_threads))$predictions)
  elapsed <- proc.time()[[3L]] - started
  rfag_stop_if(length(prediction) != nrow(validation) || any(!is.finite(prediction)), "RFAG_PREDICTION_INVALID")
  list(prediction = prediction, elapsed_seconds = as.numeric(elapsed))
}

rfag_metrics <- function(actual, predicted) {
  residual <- predicted - actual
  c(rmse = sqrt(mean(residual ^ 2)), mae = mean(abs(residual)))
}

rfag_compare_threads <- function(prepared, seed, num_trees = 500L, mtry_rule = "p_over_3", min_node_size = 20L) {
  p <- length(prepared$predictors)
  mtry <- rfag_resolve_mtry(p, mtry_rule)
  single <- rfag_fit_predict(prepared$training, prepared$validation, mtry, min_node_size, seed, num_trees, 1L)
  four <- rfag_fit_predict(prepared$training, prepared$validation, mtry, min_node_size, seed, num_trees, 4L)
  metric_single <- rfag_metrics(prepared$validation$bias_spo2_minus_sao2, single$prediction)
  metric_four <- rfag_metrics(prepared$validation$bias_spo2_minus_sao2, four$prediction)
  max_abs_difference <- max(abs(single$prediction - four$prediction))
  list(
    predictor_n = p,
    predictor_names = paste(prepared$predictors, collapse = "|"),
    resolved_mtry = mtry,
    single_elapsed_seconds = single$elapsed_seconds,
    four_elapsed_seconds = four$elapsed_seconds,
    max_abs_prediction_difference = max_abs_difference,
    rmse_single = unname(metric_single[["rmse"]]),
    rmse_four = unname(metric_four[["rmse"]]),
    mae_single = unname(metric_single[["mae"]]),
    mae_four = unname(metric_four[["mae"]]),
    exact_prediction_equal = identical(single$prediction, four$prediction),
    exact_metric_equal = identical(metric_single, metric_four)
  )
}

