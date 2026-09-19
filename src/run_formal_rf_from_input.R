args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L) || !args[[4L]] %in% c("paper", "smoke") ||
    (length(args) == 5L && args[[5L]] != "resume"))
  stop("Usage: COHORT ANALYSIS_INPUT_DIR OUTPUT_DIR paper|smoke [resume]", call. = FALSE)
cohort <- args[[1L]]
input <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
mode <- args[[4L]]
resume <- length(args) == 5L
if (!resume && (file.exists(output) || dir.exists(output))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (resume && !dir.exists(output)) stop("RESUME_OUTPUT_MISSING", call. = FALSE)
if (!cohort %in% c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang"))
  stop("COHORT_INVALID", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)
suppressPackageStartupMessages({library(digest); library(ranger); library(fastshap)})
for (path in c("code/formal/r211e_core.R", "code/cores/r211f_core.R",
               "code/formal/stage405_rf_amended_fitting_core.R",
               "code/formal/stage405_rf_amended_formal_core.R",
               "code/cores/rf_integer_oxygen_background_core.R"))
  source(file.path(root, path))
candidate_path <- file.path(input, "rf_candidate_60.rds")
receipt_path <- file.path(input, "receipt.tsv")
if (!file.exists(candidate_path) || !file.exists(receipt_path))
  stop("RF_CANDIDATE_OR_RECEIPT_MISSING", call. = FALSE)
receipt <- read.delim(receipt_path, check.names = FALSE)
row <- receipt[receipt$file == "rf_candidate_60.rds", , drop = FALSE]
if (nrow(row) != 1L || tolower(row$sha256[[1L]]) != tolower(digest::digest(file = candidate_path, algo = "sha256")))
  stop("RF_CANDIDATE_HASH_MISMATCH", call. = FALSE)
candidate <- as.data.frame(readRDS(candidate_path))
r211f_validate_candidate(candidate, cohort)
allocation <- r211e_build_allocation(candidate, cohort)
aligned <- rff_validate_allocation(candidate, allocation, cohort)
if (length(unique(candidate$analysis_patient_key)) < 30L)
  stop("RF_PATIENT_COUNT_INSUFFICIENT", call. = FALSE)
settings <- if (mode == "paper") list(trees = 500L, threads = 4L, boot = 2000L, shap_nsim = 20L) else
  list(trees = 20L, threads = 1L, boot = 20L, shap_nsim = 2L)
cohort_index <- match(cohort, r211e_expected_cohorts())
identity <- list(cohort = cohort, mode = mode, settings = settings,
                 candidate_sha256 = digest::digest(file = candidate_path, algo = "sha256"),
                 code_sha256 = vapply(c("src/run_formal_rf_from_input.R", "code/formal/r211e_core.R",
                                          "code/cores/r211f_core.R", "code/formal/stage405_rf_amended_fitting_core.R",
                                          "code/formal/stage405_rf_amended_formal_core.R",
                                          "code/cores/rf_integer_oxygen_background_core.R"),
                                        function(path) digest::digest(file = file.path(root, path), algo = "sha256"),
                                        character(1)))
identity_path <- file.path(output, "run_identity.rds")
if (resume) {
  if (!file.exists(identity_path) || !identical(readRDS(identity_path), identity))
    stop("RESUME_IDENTITY_MISMATCH", call. = FALSE)
  previous_status <- readLines(file.path(output, "status.txt"), warn = FALSE)
  if (grepl("COMPLETE", previous_status[[1L]], fixed = TRUE))
    stop("RESUME_ALREADY_COMPLETE", call. = FALSE)
} else {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output)) stop("OUTPUT_CREATE_FAILED", call. = FALSE)
if (!resume) saveRDS(identity, identity_path)
writeLines("RF_RUNNING_CANDIDATE_ONLY", file.path(output, "status.txt"))
write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", row.names = FALSE,
                                                  quote = FALSE, na = "")
tryCatch({
  source_paths <- c(candidate_path, file.path(root, "src", "run_formal_rf_from_input.R"),
                    file.path(root, "code", "formal", "r211e_core.R"),
                    file.path(root, "code", "cores", "r211f_core.R"),
                    file.path(root, "code", "formal", "stage405_rf_amended_fitting_core.R"),
                    file.path(root, "code", "formal", "stage405_rf_amended_formal_core.R"),
                    file.path(root, "code", "cores", "rf_integer_oxygen_background_core.R"))
  provenance <- data.frame(file = normalizePath(source_paths, winslash = "/"),
                           sha256 = vapply(source_paths, digest::digest, character(1),
                                           file = TRUE, algo = "sha256"))
  write_tsv(provenance, file.path(output, "input_code_hashes.tsv"))
  write_tsv(allocation, file.path(output, "allocation.tsv"))
  write_tsv(r211e_allocation_summary(candidate, allocation, cohort),
            file.path(output, "allocation_summary.tsv"))
  grid <- rff_tuning_grid()
  checkpoint_dir <- file.path(output, "checkpoints")
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir)
  fold_rows <- vector("list", 60L)
  recipe_rows <- vector("list", 6L)
  position <- 0L
  for (fold in seq_len(5L)) {
    training_index <- aligned$outer_partition == "development" & aligned$inner_fold != fold
    validation_index <- aligned$outer_partition == "development" & aligned$inner_fold == fold
    training_raw <- candidate[training_index, , drop = FALSE]
    validation_raw <- candidate[validation_index, , drop = FALSE]
    if (length(intersect(training_raw$analysis_patient_key, validation_raw$analysis_patient_key)))
      stop("RF_INNER_PATIENT_LEAKAGE", call. = FALSE)
    recipe <- r211f_fit_training_recipe(training_raw, cohort)
    r211f_validate_recipe(recipe, cohort)
    training_model <- r211f_apply_fixed_recipe(training_raw, recipe, cohort)
    validation_model <- r211f_apply_fixed_recipe(validation_raw, recipe, cohort)
    predictors <- rfag_predictors(training_model)
    if (!identical(predictors, rfag_predictors(validation_model)))
      stop("RF_INNER_PREDICTOR_MISMATCH", call. = FALSE)
    spec <- rfag_factor_spec(recipe)
    training <- rfp_transform_model_frame(rfag_prepare_frame(training_model, predictors, spec))
    validation <- rfp_transform_model_frame(rfag_prepare_frame(validation_model, predictors, spec))
    recipe_rows[[fold]] <- cbind(recipe_context = paste0("inner_fold_", fold), recipe)
    for (grid_index in seq_len(nrow(grid))) {
      position <- position + 1L
      rule <- grid$mtry_rule[[grid_index]]
      node <- as.integer(grid$min_node_size[[grid_index]])
      seed <- as.integer(20260829L + grid_index * 100L + fold)
      mtry <- rff_resolve_mtry(length(predictors), rule)
      checkpoint_path <- file.path(checkpoint_dir, sprintf("inner_f%02d_g%02d.rds", fold, grid_index))
      checkpoint_identity <- list(run = identity, fold = fold, grid_index = grid_index,
                                  training_pairs = as.character(training_raw$pair_uid),
                                  validation_pairs = as.character(validation_raw$pair_uid),
                                  recipe = recipe)
      if (file.exists(checkpoint_path)) {
        saved <- readRDS(checkpoint_path)
        if (!identical(saved$identity, checkpoint_identity)) stop("INNER_CHECKPOINT_IDENTITY_MISMATCH", call. = FALSE)
        fold_rows[[position]] <- saved$result
        next
      }
      fit <- rff_fit_model(training, mtry, node, seed, settings$trees, settings$threads)
      predicted <- rff_predict_model(fit, validation, settings$threads)
      metric <- rff_metrics(validation$bias_spo2_minus_sao2, predicted)
      one <- data.frame(
        cohort = cohort, grid_index = grid_index, inner_validation_fold = fold,
        mtry_rule = rule, min_node_size = node, predictor_n = length(predictors),
        resolved_mtry = mtry, validation_pair_n = nrow(validation),
        validation_patient_n = length(unique(validation_raw$analysis_patient_key)),
        rmse = unname(metric[["rmse"]]), mae = unname(metric[["mae"]]),
        predictor_names = paste(predictors, collapse = "|"))
      saveRDS(list(identity = checkpoint_identity, result = one), checkpoint_path)
      fold_rows[[position]] <- one
      test_interrupt <- suppressWarnings(as.integer(Sys.getenv("R9_RF_TEST_INTERRUPT_AFTER", "")))
      if (mode == "smoke" && is.finite(test_interrupt) && position == test_interrupt)
        stop("INTENTIONAL_SMOKE_INTERRUPT_AFTER_CHECKPOINT", call. = FALSE)
    }
  }
  fold_metrics <- do.call(rbind, fold_rows)
  tuning <- rff_summarize_tuning(fold_metrics, cohort)
  selected <- tuning$selected
  write_tsv(fold_metrics, file.path(output, "rf_inner_fold_metrics.tsv"))
  write_tsv(tuning$summary, file.path(output, "rf_tuning_grid_summary.tsv"))
  write_tsv(selected, file.path(output, "rf_selected_hyperparameters.tsv"))
  development_index <- aligned$outer_partition == "development"
  test_index <- aligned$outer_partition == "test"
  development_raw <- candidate[development_index, , drop = FALSE]
  test_raw <- candidate[test_index, , drop = FALSE]
  if (length(intersect(development_raw$analysis_patient_key, test_raw$analysis_patient_key)))
    stop("RF_OUTER_PATIENT_LEAKAGE", call. = FALSE)
  recipe <- r211f_fit_training_recipe(development_raw, cohort)
  r211f_validate_recipe(recipe, cohort)
  recipe_rows[[6L]] <- cbind(recipe_context = "outer_development", recipe)
  write_tsv(do.call(rbind, recipe_rows), file.path(output, "rf_training_recipes.tsv"))
  development_model <- r211f_apply_fixed_recipe(development_raw, recipe, cohort)
  test_model <- r211f_apply_fixed_recipe(test_raw, recipe, cohort)
  predictors <- rfag_predictors(development_model)
  if (!identical(predictors, rfag_predictors(test_model))) stop("RF_OUTER_PREDICTOR_MISMATCH", call. = FALSE)
  spec <- rfag_factor_spec(recipe)
  development <- rfp_transform_model_frame(rfag_prepare_frame(development_model, predictors, spec))
  test <- rfp_transform_model_frame(rfag_prepare_frame(test_model, predictors, spec))
  final_mtry <- rff_resolve_mtry(length(predictors), selected$mtry_rule[[1L]])
  final_seed <- as.integer(20260931L + cohort_index)
  final_checkpoint_path <- file.path(checkpoint_dir, "final_model_and_predictions.rds")
  final_identity <- list(run = identity, selected = selected$grid_index[[1L]],
                         development_pairs = as.character(development_raw$pair_uid),
                         test_pairs = as.character(test_raw$pair_uid), recipe = recipe)
  if (file.exists(final_checkpoint_path)) {
    saved_final <- readRDS(final_checkpoint_path)
    if (!identical(saved_final$identity, final_identity)) stop("FINAL_CHECKPOINT_IDENTITY_MISMATCH", call. = FALSE)
    fit <- saved_final$model
    predicted <- saved_final$predictions
  } else {
    fit <- rff_fit_model(development, final_mtry, selected$min_node_size[[1L]],
                         final_seed, settings$trees, settings$threads)
    predicted <- rff_predict_model(fit, test, settings$threads)
    saveRDS(list(identity = final_identity, model = fit, predictions = predicted), final_checkpoint_path)
  }
  actual <- test$bias_spo2_minus_sao2
  boot <- rff_cluster_bootstrap(actual, predicted, test_raw$analysis_patient_key,
                                 settings$boot, as.integer(20261031L + cohort_index))
  performance <- boot$summary
  performance$cohort <- cohort
  performance$development_pair_n <- nrow(development)
  performance$development_patient_n <- length(unique(development_raw$analysis_patient_key))
  performance$test_pair_n <- nrow(test)
  performance$test_patient_n <- length(unique(test_raw$analysis_patient_key))
  write_tsv(performance, file.path(output, "rf_outer_test_performance.tsv"))
  write_tsv(data.frame(pair_uid = test_raw$pair_uid, analysis_patient_key = test_raw$analysis_patient_key,
                       observed_bias = actual, predicted_bias = predicted),
            file.path(output, "rf_outer_test_predictions.tsv"))
  saveRDS(fit, file.path(output, "rf_final_model.rds"))
  saveRDS(boot$draws, file.path(output, "rf_patient_bootstrap_draws.rds"))
  shap_index <- rff_patient_balanced_shap_indices(test_raw$analysis_patient_key,
                                                    test_raw$pair_uid, 1000L,
                                                    as.integer(20261131L + cohort_index))
  shap_x <- test[shap_index, predictors, drop = FALSE]
  set.seed(as.integer(20261231L + cohort_index))
  shap_matrix <- as.matrix(fastshap::explain(fit, X = shap_x,
      pred_wrapper = function(object, newdata) rff_predict_model(object, newdata, settings$threads),
      nsim = settings$shap_nsim, adjust = TRUE, parallel = FALSE))
  if (!identical(colnames(shap_matrix), predictors) || any(!is.finite(shap_matrix)))
    stop("RF_SHAP_INVALID", call. = FALSE)
  grouped <- rff_group_shap(shap_matrix)
  write_tsv(grouped$summary, file.path(output, "rf_shap_summary.tsv"))
  write_tsv(grouped$model_column_map, file.path(output, "rf_shap_model_column_map.tsv"))
  write_tsv(data.frame(pair_uid = test_raw$pair_uid[shap_index],
                       analysis_patient_key = test_raw$analysis_patient_key[shap_index]),
            file.path(output, "rf_shap_sample.tsv"))
  display_features <- c("pair_uid", "age_value_per_10y", "age_interval", "sex_class",
                        "cci_value", "sofa_total", "pao2_candidate_value",
                        "paco2_candidate_value", "ph_candidate_value",
                        "lactate_candidate_value", "hb_candidate_value",
                        "concurrent_vasoactive_medication_use_state",
                        "current_invasive_mechanical_ventilation_state")
  essential_display_features <- c("pair_uid", "age_interval", "sex_class",
                                  "pao2_candidate_value", "paco2_candidate_value",
                                  "ph_candidate_value")
  if (!all(essential_display_features %in% names(test_raw)) ||
      !"paired_mean_saturation_c90__model" %in% names(test)) {
    stop("RF_SHAP_DISPLAY_FEATURES_MISSING", call. = FALSE)
  }
  feature_values <- test_raw[shap_index, intersect(display_features, names(test_raw)), drop = FALSE]
  feature_values$paired_mean_saturation_c90_rf <-
    test$paired_mean_saturation_c90__model[shap_index]
  write_tsv(feature_values, file.path(output, "rf_shap_feature_values.tsv"))
  saveRDS(shap_matrix, file.path(output, "rf_shap_model_columns.rds"))
  saveRDS(grouped$values, file.path(output, "rf_shap_grouped.rds"))
  write_tsv(data.frame(cohort = cohort, mode = mode, trees = settings$trees,
                       threads = settings$threads, bootstrap_replicates = settings$boot,
                       shap_permutations = settings$shap_nsim, selected_grid = selected$grid_index,
                       predictor_n = length(predictors), predictors = paste(predictors, collapse = "|"),
                       paired_mean_integerized = all(development$paired_mean_saturation_c90__model %% 1 == 0),
                       standalone_spo2_or_sao2 = any(grepl("^(spo2|sao2)_percent", predictors))),
            file.path(output, "rf_diagnostics.tsv"))
  writeLines(if (mode == "paper") "RF_PAPER_SETTINGS_COMPLETE_CANDIDATE_ONLY" else
               "RF_SMOKE_COMPLETE_NOT_PAPER_SETTINGS", file.path(output, "status.txt"))
  cat("RF_SYNTHETIC_OR_FORMAL_INPUT_RUN_PASS mode=", mode, "\n")
}, error = function(e) {
  writeLines(paste("RF_FAIL", conditionMessage(e)), file.path(output, "status.txt"))
  stop(e)
})
