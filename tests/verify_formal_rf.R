args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: ANALYSIS_INPUT_DIR RF_OUTPUT_DIR NEW_QA_DIR", call. = FALSE)
input <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
if (file.exists(qa) || dir.exists(qa)) stop("QA_OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
checks <- character()
check <- function(condition, label) {
  if (!isTRUE(condition)) stop(paste0("RF_INDEPENDENT_QA_FAIL:", label), call. = FALSE)
  checks <<- c(checks, label)
}
read_tsv <- function(name) read.delim(file.path(result, name), check.names = FALSE,
                                     stringsAsFactors = FALSE)
tryCatch({
  library(digest)
  candidate_path <- file.path(input, "rf_candidate_60.rds")
  candidate <- as.data.frame(readRDS(candidate_path))
  hashes <- read_tsv("input_code_hashes.tsv")
  check(all(file.exists(hashes$file)), "all_declared_inputs_exist")
  check(all(tolower(vapply(hashes$file, digest::digest, character(1),
                           file = TRUE, algo = "sha256")) == tolower(hashes$sha256)),
        "declared_file_hashes_match")
  allocation <- read_tsv("allocation.tsv")
  check(nrow(allocation) == nrow(candidate), "allocation_covers_all_pairs")
  check(!anyDuplicated(allocation$pair_uid), "allocation_pair_unique")
  check(setequal(allocation$pair_uid, candidate$pair_uid), "allocation_pair_ids_match")
  map <- unique(allocation[, c("analysis_patient_key", "outer_partition")])
  check(!anyDuplicated(map$analysis_patient_key), "one_outer_partition_per_patient")
  check(setequal(map$outer_partition, c("development", "test")), "both_outer_partitions")
  check(all(is.na(allocation$inner_fold[allocation$outer_partition == "test"])),
        "test_has_no_inner_fold")
  check(setequal(na.omit(unique(allocation$inner_fold)), seq_len(5L)), "five_inner_folds")
  fold_map <- unique(allocation[allocation$outer_partition == "development",
                                c("analysis_patient_key", "inner_fold")])
  check(!anyDuplicated(fold_map$analysis_patient_key), "one_fold_per_development_patient")
  check(!any(grepl("^(spo2|sao2)_percent__(model|was_missing|was_unseen)$",
                   names(candidate))), "no_oxygen_model_fields_in_candidate")
  check(isTRUE(all.equal(candidate$bias_spo2_minus_sao2,
                         candidate$spo2_percent - candidate$sao2_percent)), "bias_identity")
  check(isTRUE(all.equal(candidate$paired_mean_saturation_c90,
                         (candidate$spo2_percent + candidate$sao2_percent) / 2 - 90)),
        "paired_mean_identity")
  recipes <- read_tsv("rf_training_recipes.tsv")
  check(setequal(recipes$recipe_context, c(paste0("inner_fold_", 1:5), "outer_development")),
        "six_training_only_recipes")
  for (context in unique(recipes$recipe_context)) {
    x <- recipes[recipes$recipe_context == context, , drop = FALSE]
    expected <- if (context == "outer_development")
      sum(allocation$outer_partition == "development") else {
        fold <- as.integer(sub("inner_fold_", "", context))
        sum(allocation$outer_partition == "development" & allocation$inner_fold != fold)
      }
    structural <- x$feature_role == "structural"
    check(all(x$observed_n[structural] == expected & x$retained[structural]),
          paste0(context, "_structural_always_retained"))
    clinical <- x$feature_role == "candidate" & x$retained
    check(all(x$observed_n[clinical] >= 500 & x$observed_support[clinical] >= 0.10),
          paste0(context, "_clinical_support_threshold"))
    check(!any(x$feature %in% c("spo2_percent", "sao2_percent", "bias_spo2_minus_sao2")),
          paste0(context, "_no_leakage_features"))
  }
  folds <- read_tsv("rf_inner_fold_metrics.tsv")
  check(nrow(folds) == 60L, "sixty_inner_tasks")
  check(nrow(unique(folds[, c("grid_index", "inner_validation_fold")])) == 60L,
        "inner_task_identity_unique")
  check(all(is.finite(folds$rmse) & is.finite(folds$mae)), "inner_metrics_finite")
  selected <- read_tsv("rf_selected_hyperparameters.tsv")
  summary <- read_tsv("rf_tuning_grid_summary.tsv")
  check(nrow(selected) == 1L & nrow(summary) == 12L, "tuning_selection_present")
  ranked <- summary[order(summary$mean_rmse, summary$mean_mae,
                          -summary$min_node_size,
                          match(summary$mtry_rule, c("p_over_4", "p_over_3", "p_over_2", "p"))), ]
  check(selected$grid_index[[1L]] == ranked$grid_index[[1L]], "selected_rule_recomputed")
  diagnostics <- read_tsv("rf_diagnostics.tsv")
  check(nrow(diagnostics) == 1L & !diagnostics$standalone_spo2_or_sao2[[1L]],
        "no_standalone_oxygen_predictors")
  check(diagnostics$paired_mean_integerized[[1L]], "paired_mean_integerized")
  checkpoint_paths <- list.files(file.path(result, "checkpoints"),
                                 pattern = "^inner_f[0-9]{2}_g[0-9]{2}\\.rds$", full.names = TRUE)
  check(length(checkpoint_paths) == 60L, "sixty_resumable_inner_checkpoints")
  run_identity <- readRDS(file.path(result, "run_identity.rds"))
  check(all(vapply(checkpoint_paths, function(path)
    identical(readRDS(path)$identity$run, run_identity), logical(1))),
    "checkpoint_run_identity_matches")
  final_checkpoint <- readRDS(file.path(result, "checkpoints", "final_model_and_predictions.rds"))
  check(identical(final_checkpoint$identity$run, run_identity),
        "final_checkpoint_run_identity_matches")
  if (diagnostics$mode[[1L]] == "paper")
    check(diagnostics$trees[[1L]] == 500L & diagnostics$bootstrap_replicates[[1L]] == 2000L &
            diagnostics$shap_permutations[[1L]] == 20L, "paper_configuration")
  predictions <- read_tsv("rf_outer_test_predictions.tsv")
  check(nrow(predictions) == sum(allocation$outer_partition == "test"), "test_predictions_complete")
  check(!anyDuplicated(predictions$pair_uid), "test_predictions_unique")
  check(all(is.finite(predictions$observed_bias) & is.finite(predictions$predicted_bias)),
        "test_predictions_finite")
  test_raw <- candidate[match(predictions$pair_uid, candidate$pair_uid), ]
  check(!anyNA(test_raw$pair_uid), "test_pairs_found_in_candidate")
  check(isTRUE(all.equal(predictions$observed_bias, test_raw$bias_spo2_minus_sao2)),
        "test_observed_bias_matches_input")
  perf <- read_tsv("rf_outer_test_performance.tsv")
  delta <- predictions$predicted_bias - predictions$observed_bias
  rmse <- sqrt(mean(delta^2))
  mae <- mean(abs(delta))
  check(abs(perf$estimate[perf$metric == "rmse"] - rmse) < 1e-10,
        "test_rmse_recomputed")
  check(abs(perf$estimate[perf$metric == "mae"] - mae) < 1e-10,
        "test_mae_recomputed")
  draws <- readRDS(file.path(result, "rf_patient_bootstrap_draws.rds"))
  check(nrow(draws) == diagnostics$bootstrap_replicates[[1L]], "bootstrap_draw_count")
  for (metric in c("rmse", "mae")) {
    x <- perf[perf$metric == metric, ]
    bounds <- quantile(draws[, metric], probs = c(0.025, 0.975), na.rm = TRUE)
    check(max(abs(c(x$ci_lower, x$ci_upper) - as.numeric(bounds))) < 1e-10,
          paste0(metric, "_bootstrap_ci_recomputed"))
  }
  sample <- read_tsv("rf_shap_sample.tsv")
  check(!anyDuplicated(sample$analysis_patient_key), "shap_one_pair_per_patient")
  check(all(sample$pair_uid %in% predictions$pair_uid), "shap_sample_from_test_only")
  model_shap <- as.matrix(readRDS(file.path(result, "rf_shap_model_columns.rds")))
  grouped_shap <- as.matrix(readRDS(file.path(result, "rf_shap_grouped.rds")))
  check(nrow(model_shap) == nrow(sample), "shap_rows_match_sample")
  check(max(abs(rowSums(model_shap) - rowSums(grouped_shap))) < 1e-8,
        "shap_grouping_additive")
  shap_summary <- read_tsv("rf_shap_summary.tsv")
  group_index <- match(shap_summary$feature_group, colnames(grouped_shap))
  check(!anyNA(group_index), "shap_groups_present")
  check(max(abs(shap_summary$mean_abs_shap - colMeans(abs(grouped_shap)) [group_index])) < 1e-10,
        "shap_importance_recomputed")
  writeLines(c("RF_INDEPENDENT_QA_PASS", paste0("checks=", length(checks)), checks),
             file.path(qa, "qa_status.txt"))
  cat("RF_INDEPENDENT_QA_PASS checks=", length(checks), "\n")
}, error = function(e) {
  writeLines(paste("RF_INDEPENDENT_QA_FAIL", conditionMessage(e)),
             file.path(qa, "qa_status.txt"))
  stop(e)
})
