rff_stop_if <- function(condition, code) {
  if (isTRUE(condition)) stop(code, call. = FALSE)
  invisible(TRUE)
}

rff_hash_file <- function(path) {
  rff_stop_if(!file.exists(path), paste0("RFF_FILE_MISSING:", path))
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

rff_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid(), "-", sprintf("%08d", sample.int(1e8, 1L)))
  saveRDS(object, temporary, compress = FALSE)
  moved <- file.rename(temporary, path)
  if (!moved) {
    copied <- file.copy(temporary, path, overwrite = TRUE)
    unlink(temporary)
    rff_stop_if(!copied, "RFF_ATOMIC_RDS_REPLACE_FAILED")
  }
  invisible(path)
}

rff_atomic_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid(), "-", sprintf("%08d", sample.int(1e8, 1L)))
  utils::write.table(x, temporary, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  moved <- file.rename(temporary, path)
  if (!moved) {
    copied <- file.copy(temporary, path, overwrite = TRUE)
    unlink(temporary)
    rff_stop_if(!copied, "RFF_ATOMIC_TSV_REPLACE_FAILED")
  }
  invisible(path)
}

rff_atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid(), "-", sprintf("%08d", sample.int(1e8, 1L)))
  writeLines(as.character(x), temporary, useBytes = TRUE)
  moved <- file.rename(temporary, path)
  if (!moved) {
    copied <- file.copy(temporary, path, overwrite = TRUE)
    unlink(temporary)
    rff_stop_if(!copied, "RFF_ATOMIC_TEXT_REPLACE_FAILED")
  }
  invisible(path)
}

rff_atomic_write_parquet <- function(x, path) {
  rff_stop_if(!requireNamespace("arrow", quietly = TRUE), "RFF_ARROW_UNAVAILABLE")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid(), "-", sprintf("%08d", sample.int(1e8, 1L)))
  arrow::write_parquet(x, temporary, compression = "zstd")
  moved <- file.rename(temporary, path)
  if (!moved) {
    copied <- file.copy(temporary, path, overwrite = TRUE)
    unlink(temporary)
    rff_stop_if(!copied, "RFF_ATOMIC_PARQUET_REPLACE_FAILED")
  }
  invisible(path)
}

rff_validate_output_root <- function(project_root, output_root) {
  project <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  output <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  approved <- normalizePath(file.path(project, "07_analysis", "05_runs"), winslash = "/", mustWork = FALSE)
  rff_stop_if(!startsWith(tolower(output), paste0(tolower(approved), "/")), "RFF_OUTPUT_OUTSIDE_APPROVED_RUN_ROOT")
  rff_stop_if(grepl("\\$", output), "RFF_OUTPUT_UNRESOLVED_VARIABLE")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  normalizePath(output, winslash = "/", mustWork = TRUE)
}

rff_assert_file_hash <- function(path, expected_hash, code = "RFF_INPUT_HASH_MISMATCH") {
  rff_stop_if(!identical(tolower(rff_hash_file(path)), tolower(as.character(expected_hash))), code)
  invisible(TRUE)
}

rff_fit_model <- function(training, mtry, min_node_size, seed, num_trees = 500L, num_threads = 4L) {
  rff_stop_if(!requireNamespace("ranger", quietly = TRUE), "RFF_RANGER_UNAVAILABLE")
  ranger::ranger(
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
}

rff_predict_model <- function(model, newdata, num_threads = 4L) {
  prediction <- as.numeric(stats::predict(model, data = newdata, num.threads = as.integer(num_threads))$predictions)
  rff_stop_if(length(prediction) != nrow(newdata) || any(!is.finite(prediction)), "RFF_PREDICTION_INVALID")
  prediction
}

rff_mtry_rules <- function() c("p_over_4", "p_over_3", "p_over_2", "p")

rff_tuning_grid <- function() {
  expand.grid(
    mtry_rule = rff_mtry_rules(),
    min_node_size = c(5L, 20L, 50L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

rff_resolve_mtry <- function(p, rule) {
  divisors <- c(p_over_4 = 4, p_over_3 = 3, p_over_2 = 2, p = 1)
  rff_stop_if(length(p) != 1L || !is.finite(p) || p < 1L || !rule %in% names(divisors), "RFF_MTRY_INPUT_INVALID")
  as.integer(max(1L, min(as.integer(p), round(as.integer(p) / divisors[[rule]]))))
}

rff_summarize_tuning <- function(fold_metrics, cohort) {
  required <- c("cohort", "grid_index", "inner_validation_fold", "mtry_rule", "min_node_size", "predictor_n", "resolved_mtry", "rmse", "mae")
  rff_stop_if(!all(required %in% names(fold_metrics)), "RFF_TUNING_INTERFACE_MISSING")
  expected <- expand.grid(grid_index = seq_len(12L), inner_validation_fold = seq_len(5L))
  observed <- unique(fold_metrics[, c("grid_index", "inner_validation_fold")])
  rff_stop_if(nrow(fold_metrics) != 60L || nrow(observed) != 60L || !setequal(do.call(paste, expected), do.call(paste, observed)), "RFF_TUNING_TASK_CLOSURE_INVALID")
  rff_stop_if(!all(fold_metrics$cohort == cohort), "RFF_TUNING_COHORT_INVALID")
  rff_stop_if(any(!is.finite(fold_metrics$rmse)) || any(!is.finite(fold_metrics$mae)), "RFF_TUNING_METRIC_INVALID")
  grid <- rff_tuning_grid()
  summary_rows <- lapply(seq_len(12L), function(index) {
    x <- fold_metrics[fold_metrics$grid_index == index, , drop = FALSE]
    rff_stop_if(nrow(x) != 5L || length(unique(x$mtry_rule)) != 1L || length(unique(x$min_node_size)) != 1L, "RFF_GRID_METADATA_INVALID")
    data.frame(
      cohort = cohort,
      grid_index = index,
      mtry_rule = unique(x$mtry_rule),
      min_node_size = unique(x$min_node_size),
      predictor_n_min = min(x$predictor_n),
      predictor_n_max = max(x$predictor_n),
      resolved_mtry_min = min(x$resolved_mtry),
      resolved_mtry_max = max(x$resolved_mtry),
      mean_rmse = mean(x$rmse),
      mean_mae = mean(x$mae),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary_rows)
  rff_stop_if(!identical(as.character(summary$mtry_rule), as.character(grid$mtry_rule)) || !identical(as.integer(summary$min_node_size), as.integer(grid$min_node_size)), "RFF_GRID_ORDER_INVALID")
  rule_rank <- match(summary$mtry_rule, rff_mtry_rules())
  order_index <- order(summary$mean_rmse, summary$mean_mae, -summary$min_node_size, rule_rank)
  list(summary = summary, selected = summary[order_index[[1L]], , drop = FALSE])
}

rff_metrics <- function(actual, predicted) {
  rff_stop_if(length(actual) != length(predicted) || length(actual) < 2L || any(!is.finite(actual)) || any(!is.finite(predicted)), "RFF_METRIC_INPUT_INVALID")
  residual <- predicted - actual
  denominator <- sum((actual - mean(actual)) ^ 2)
  prediction_variance <- stats::var(predicted)
  calibration <- if (is.finite(prediction_variance) && prediction_variance > 0) {
    unname(stats::coef(stats::lm(actual ~ predicted)))
  } else {
    c(NA_real_, NA_real_)
  }
  c(
    rmse = sqrt(mean(residual ^ 2)),
    mae = mean(abs(residual)),
    r_squared = if (denominator > 0) 1 - sum(residual ^ 2) / denominator else NA_real_,
    calibration_intercept = calibration[[1L]],
    calibration_slope = calibration[[2L]]
  )
}

rff_cluster_bootstrap <- function(actual, predicted, patient, b = 2000L, seed) {
  rff_stop_if(length(actual) != length(predicted) || length(actual) != length(patient), "RFF_BOOTSTRAP_LENGTH_INVALID")
  rff_stop_if(anyNA(patient) || length(unique(patient)) < 2L || b < 1L, "RFF_BOOTSTRAP_CLUSTER_INVALID")
  set.seed(as.integer(seed))
  patient <- as.character(patient)
  by_patient <- split(seq_along(patient), patient)
  ids <- sort(names(by_patient))
  metric_names <- names(rff_metrics(actual, predicted))
  draws <- matrix(NA_real_, nrow = b, ncol = length(metric_names), dimnames = list(NULL, metric_names))
  for (iteration in seq_len(b)) {
    sampled <- sample(ids, length(ids), replace = TRUE)
    index <- unlist(by_patient[sampled], use.names = FALSE)
    draws[iteration, ] <- rff_metrics(actual[index], predicted[index])
  }
  successful <- colSums(is.finite(draws))
  list(
    summary = data.frame(
      metric = colnames(draws),
      estimate = unname(rff_metrics(actual, predicted)),
      ci_lower = apply(draws, 2L, stats::quantile, probs = 0.025, na.rm = TRUE),
      ci_upper = apply(draws, 2L, stats::quantile, probs = 0.975, na.rm = TRUE),
      bootstrap_requested = as.integer(b),
      bootstrap_successful = as.integer(successful),
      stringsAsFactors = FALSE
    ),
    draws = draws
  )
}

rff_patient_balanced_shap_indices <- function(patient, pair_uid, maximum_patients = 1000L, seed) {
  rff_stop_if(length(patient) != length(pair_uid) || length(patient) == 0L || anyNA(patient) || anyNA(pair_uid), "RFF_SHAP_SAMPLE_INPUT_INVALID")
  patient <- as.character(patient)
  pair_uid <- as.character(pair_uid)
  rff_stop_if(anyDuplicated(pair_uid) > 0L, "RFF_SHAP_PAIR_DUPLICATED")
  ids <- sort(unique(patient))
  set.seed(as.integer(seed))
  selected_patients <- if (length(ids) <= maximum_patients) ids else sample(ids, maximum_patients, replace = FALSE)
  selected_patients <- sort(selected_patients)
  indices <- vapply(selected_patients, function(id) {
    candidates <- which(patient == id)
    candidates <- candidates[order(pair_uid[candidates])]
    candidates[[sample.int(length(candidates), 1L)]]
  }, integer(1L))
  rff_stop_if(length(indices) != min(length(ids), maximum_patients) || anyDuplicated(patient[indices]) > 0L, "RFF_SHAP_SAMPLE_BALANCE_INVALID")
  indices
}

rff_feature_group <- function(model_column) {
  base <- sub("__(model|was_missing|was_unseen)$", "", model_column)
  if (base %in% c("age_value_per_10y", "age_topcoded_indicator", "age_interval")) return("Age")
  labels <- c(
    paired_mean_saturation_c90 = "Paired mean saturation",
    sex_class = "Sex",
    cci_value = "Charlson Comorbidity Index",
    sofa_total = "SOFA",
    pao2_candidate_value = "PaO2",
    paco2_candidate_value = "PaCO2",
    ph_candidate_value = "pH",
    lactate_candidate_value = "Lactate",
    hb_candidate_value = "Hemoglobin",
    concurrent_vasoactive_medication_use_state = "Concurrent vasoactive medication use",
    current_invasive_mechanical_ventilation_state = "Current invasive mechanical ventilation"
  )
  if (base %in% names(labels)) unname(labels[[base]]) else base
}

rff_group_shap <- function(shap_matrix) {
  shap_matrix <- as.matrix(shap_matrix)
  rff_stop_if(is.null(colnames(shap_matrix)) || any(!is.finite(shap_matrix)), "RFF_SHAP_MATRIX_INVALID")
  groups <- vapply(colnames(shap_matrix), rff_feature_group, character(1L))
  group_order <- unique(groups)
  grouped <- sapply(group_order, function(group) rowSums(shap_matrix[, groups == group, drop = FALSE]))
  if (is.null(dim(grouped))) grouped <- matrix(grouped, ncol = 1L, dimnames = list(NULL, group_order))
  colnames(grouped) <- group_order
  rff_stop_if(max(abs(rowSums(grouped) - rowSums(shap_matrix))) > 1e-10, "RFF_GROUPED_SHAP_NONADDITIVE")
  list(
    values = grouped,
    summary = data.frame(
      feature_group = colnames(grouped),
      feature_role = ifelse(colnames(grouped) == "Paired mean saturation", "structural_oxygenation_background", "clinical_characteristic"),
      mean_abs_shap = colMeans(abs(grouped)),
      mean_shap = colMeans(grouped),
      stringsAsFactors = FALSE
    ),
    model_column_map = data.frame(model_column = colnames(shap_matrix), feature_group = groups, stringsAsFactors = FALSE)
  )
}

rff_checkpoint_identity <- function(cohort, fold, grid_index, seed, input_hashes, code_hashes, num_trees = 500L, num_threads = 4L) {
  list(
    cohort = as.character(cohort),
    fold = as.integer(fold),
    grid_index = as.integer(grid_index),
    seed = as.integer(seed),
    input_hashes = as.character(input_hashes),
    code_hashes = as.character(code_hashes),
    num_trees = as.integer(num_trees),
    num_threads = as.integer(num_threads)
  )
}

rff_checkpoint_valid <- function(checkpoint, expected_identity) {
  is.list(checkpoint) && !is.null(checkpoint$identity) && !is.null(checkpoint$result) && identical(checkpoint$identity, expected_identity)
}

rff_outer_access_token <- function(fold_metrics, selected, cohort, input_hashes, code_hashes) {
  summary <- rff_summarize_tuning(fold_metrics, cohort)
  rff_stop_if(nrow(selected) != 1L || !identical(as.integer(selected$grid_index), as.integer(summary$selected$grid_index)), "RFF_OUTER_SELECTED_RULE_NOT_FROZEN")
  selected_key <- data.frame(
    grid_index = as.integer(selected$grid_index[[1L]]),
    mtry_rule = as.character(selected$mtry_rule[[1L]]),
    min_node_size = as.integer(selected$min_node_size[[1L]]),
    stringsAsFactors = FALSE
  )
  digest::digest(list(
    cohort = cohort,
    task_count = nrow(fold_metrics),
    selected = selected_key,
    input_hashes = sort(as.character(input_hashes)),
    code_hashes = sort(as.character(code_hashes))
  ), algo = "sha256", serialize = TRUE)
}

rff_validate_allocation <- function(candidate, allocation, cohort) {
  required_candidate <- c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid", "bias_spo2_minus_sao2")
  required_allocation <- c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid", "outer_partition", "inner_fold")
  rff_stop_if(!all(required_candidate %in% names(candidate)) || !all(required_allocation %in% names(allocation)), "RFF_OUTER_INTERFACE_MISSING")
  rff_stop_if(nrow(candidate) == 0L || nrow(candidate) != nrow(allocation), "RFF_OUTER_ROW_COUNT_INVALID")
  rff_stop_if(anyDuplicated(candidate$pair_uid) || anyDuplicated(allocation$pair_uid), "RFF_OUTER_PAIR_DUPLICATED")
  index <- match(as.character(candidate$pair_uid), as.character(allocation$pair_uid))
  rff_stop_if(anyNA(index), "RFF_OUTER_ALLOCATION_INCOMPLETE")
  aligned <- allocation[index, , drop = FALSE]
  rff_stop_if(!all(as.character(candidate$cohort) == cohort) || !all(as.character(aligned$cohort) == cohort), "RFF_OUTER_COHORT_INVALID")
  rff_stop_if(!identical(as.character(candidate$analysis_patient_key), as.character(aligned$analysis_patient_key)), "RFF_OUTER_PATIENT_KEY_MISMATCH")
  rff_stop_if(!identical(as.character(candidate$analysis_stay_key), as.character(aligned$analysis_stay_key)), "RFF_OUTER_STAY_KEY_MISMATCH")
  rff_stop_if(!all(as.character(aligned$outer_partition) %in% c("development", "test")), "RFF_OUTER_PARTITION_INVALID")
  development_patients <- unique(as.character(aligned$analysis_patient_key[aligned$outer_partition == "development"]))
  test_patients <- unique(as.character(aligned$analysis_patient_key[aligned$outer_partition == "test"]))
  rff_stop_if(length(intersect(development_patients, test_patients)) > 0L, "RFF_OUTER_PATIENT_LEAKAGE")
  aligned
}

rff_restore_recipe_attributes <- function(recipe, cohort) {
  attr(recipe, "cohort") <- cohort
  attr(recipe, "minimum_support") <- 0.10
  attr(recipe, "minimum_observed_n") <- 500L
  recipe
}

rff_prepare_outer <- function(candidate, allocation, recipe_bundle, cohort, outer_token, expected_outer_token) {
  rff_stop_if(!identical(as.character(outer_token), as.character(expected_outer_token)), "RFF_OUTER_ACCESS_GUARD_FAILED")
  aligned <- rff_validate_allocation(candidate, allocation, cohort)
  recipe <- recipe_bundle[as.character(recipe_bundle$recipe_context) == "outer_development", , drop = FALSE]
  rff_stop_if(nrow(recipe) == 0L || anyDuplicated(recipe$feature), "RFF_OUTER_RECIPE_INVALID")
  recipe <- rff_restore_recipe_attributes(recipe, cohort)
  development_index <- as.character(aligned$outer_partition) == "development"
  test_index <- as.character(aligned$outer_partition) == "test"
  rff_stop_if(!any(development_index) || !any(test_index), "RFF_OUTER_PARTITION_EMPTY")
  development_raw <- candidate[development_index, , drop = FALSE]
  test_raw <- candidate[test_index, , drop = FALSE]
  development_model <- r211f_apply_fixed_recipe(development_raw, recipe, cohort)
  test_model <- r211f_apply_fixed_recipe(test_raw, recipe, cohort)
  predictors_development <- rfag_predictors(development_model)
  predictors_test <- rfag_predictors(test_model)
  rff_stop_if(!identical(predictors_development, predictors_test), "RFF_OUTER_PREDICTOR_ORDER_MISMATCH")
  spec <- rfag_factor_spec(recipe)
  development_frame <- rfag_prepare_frame(development_model, predictors_development, spec)
  test_frame <- rfag_prepare_frame(test_model, predictors_development, spec)
  rff_stop_if(length(intersect(as.character(development_raw$analysis_patient_key), as.character(test_raw$analysis_patient_key))) > 0L, "RFF_OUTER_PATIENT_LEAKAGE_AFTER_PREPARATION")
  list(
    predictors = predictors_development,
    recipe = recipe,
    development = development_frame,
    test = test_frame,
    development_meta = development_raw[, c("analysis_patient_key", "analysis_stay_key", "pair_uid"), drop = FALSE],
    test_meta = test_raw[, c("analysis_patient_key", "analysis_stay_key", "pair_uid"), drop = FALSE],
    test_source_features = test_raw[, as.character(recipe$feature[recipe$retained]), drop = FALSE]
  )
}

rff_inner_checkpoint_valid <- function(path, expected_identity) {
  if (!file.exists(path)) return(list(valid = FALSE, reason = "absent", checkpoint = NULL))
  checkpoint <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(checkpoint)) return(list(valid = FALSE, reason = "corrupt", checkpoint = NULL))
  if (!rff_checkpoint_valid(checkpoint, expected_identity)) return(list(valid = FALSE, reason = "identity_mismatch", checkpoint = checkpoint))
  required <- c("cohort", "grid_index", "inner_validation_fold", "mtry_rule", "min_node_size", "predictor_n", "resolved_mtry", "rmse", "mae")
  if (!all(required %in% names(checkpoint$result))) return(list(valid = FALSE, reason = "result_interface_invalid", checkpoint = checkpoint))
  if (any(!is.finite(c(checkpoint$result$rmse, checkpoint$result$mae)))) return(list(valid = FALSE, reason = "result_value_invalid", checkpoint = checkpoint))
  list(valid = TRUE, reason = "valid", checkpoint = checkpoint)
}

rff_quarantine_file <- function(path, quarantine_root, reason) {
  if (!file.exists(path)) return(NA_character_)
  dir.create(quarantine_root, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(quarantine_root, paste0(basename(path), ".", reason, ".", format(Sys.time(), "%Y%m%dT%H%M%S")))
  rff_stop_if(!file.rename(path, target), "RFF_QUARANTINE_FAILED")
  normalizePath(target, winslash = "/", mustWork = TRUE)
}
