args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("USAGE_UPSTREAM_AND_OUTPUT_DIRECTORY", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
candidate_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
run_parent <- normalizePath(file.path(candidate_root, "runs"), winslash = "/", mustWork = TRUE)
upstream_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_root <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (!identical(dirname(upstream_root), run_parent) || !identical(dirname(output_root), run_parent)) stop("RUN_PATH_OUTSIDE_CANDIDATE", call. = FALSE)
if (dir.exists(output_root) || file.exists(output_root)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (!file.exists(file.path(upstream_root, "status.txt")) ||
    !file.exists(file.path(upstream_root, "independent_qa_status.txt")) ||
    !identical(trimws(readLines(file.path(upstream_root, "status.txt"), warn = FALSE)), "SYNTHETIC_UPSTREAM_PASS") ||
    !identical(trimws(readLines(file.path(upstream_root, "independent_qa_status.txt"), warn = FALSE)), "SYNTHETIC_UPSTREAM_INDEPENDENT_QA_PASS")) stop("UPSTREAM_GATE_NOT_PASSED", call. = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
writeLines("SYNTHETIC_RF_RUNNING", file.path(output_root, "status.txt"))
options(error = function() {
  writeLines("SYNTHETIC_RF_FAIL", file.path(output_root, "status.txt"))
  quit(save = "no", status = 1L)
})

suppressPackageStartupMessages({library(digest); library(ranger); library(fastshap)})
module_root <- file.path(candidate_root, "code", "cores")
spec <- data.frame(module = c("R211B", "R211F", "RF_FITTING_GATE", "RF_FORMAL", "RF_INTEGER_OXYGEN"), relative_path = c(
  "r211b_core.R",
  "r211f_core.R",
  "stage405_rf_amended_fitting_core.R",
  "rff_demo_core.R",
  "rf_integer_oxygen_background_core.R"
), expected_sha256 = c(
  "68F51663252389676FBC58BAF7CE34B6671F155D9A7801C66B47A829EEE99C09",
  "21B8B595BDBAD0CA6AF854771E6D70C0FDBFBDE2DBBD1884855A1D6498630E87",
  "53F5AB508FC806DF6D8F0770255FEA91DAE05839B7ADD1F8CF81D69183A23163",
  "D1E883BD47EC3FAACC9A8B51EE55D9427029F9EE39A9DDA7470172AB2A97452C",
  "D141F36A3ECBE97574877D916A207146593EC9B4829686C9290A0B1651A9EE2E"
), stringsAsFactors = FALSE)
paths <- file.path(module_root, spec$relative_path)
if (any(!file.exists(paths))) stop("MODULE_SOURCE_MISSING", call. = FALSE)
spec$observed_sha256 <- toupper(unname(tools::sha256sum(paths)))
if (any(spec$observed_sha256 != spec$expected_sha256)) stop("MODULE_SOURCE_HASH_MISMATCH", call. = FALSE)
for (path in paths) source(path, local = FALSE)
source(file.path(candidate_root, "src", "synthetic_model_inputs.R"), local = FALSE)

upstream <- readRDS(file.path(upstream_root, "synthetic_pairs_INTERNAL_ONLY.rds"))
if (!is.list(upstream) || !all(c("M60", "M5") %in% names(upstream))) stop("SYNTHETIC_PAIR_INTERFACE_INVALID", call. = FALSE)
candidate <- private_build_model_inputs(upstream$M60)$candidate
patient <- unique(as.character(candidate$analysis_patient_key))
if (length(patient) < 200L) stop("RF_FIXTURE_REQUIRES_200_SYNTHETIC_PATIENTS", call. = FALSE)
candidate <- candidate[as.character(candidate$analysis_patient_key) %in% patient[seq_len(200L)], , drop = FALSE]
patient <- patient[seq_len(200L)]
allocation <- candidate[, c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid"), drop = FALSE]
development_patient <- patient[seq_len(160L)]
allocation$outer_partition <- ifelse(allocation$analysis_patient_key %in% development_patient, "development", "test")
fold_map <- stats::setNames(rep(seq_len(5L), length.out = 160L), development_patient)
allocation$inner_fold <- as.integer(ifelse(allocation$outer_partition == "development", fold_map[allocation$analysis_patient_key], NA_integer_))
if (length(intersect(unique(allocation$analysis_patient_key[allocation$outer_partition == "development"]),
                     unique(allocation$analysis_patient_key[allocation$outer_partition == "test"]))) != 0L) stop("PATIENT_SPLIT_LEAKAGE", call. = FALSE)
if (nrow(candidate) != 800L || sum(allocation$outer_partition == "development") != 640L ||
    sum(allocation$outer_partition == "test") != 160L) stop("RF_SYNTHETIC_SPLIT_CARDINALITY", call. = FALSE)

grid <- rff_tuning_grid()
fold_rows <- vector("list", nrow(grid) * 5L)
k <- 0L
for (fold in seq_len(5L)) {
  training_index <- allocation$outer_partition == "development" & allocation$inner_fold != fold
  validation_index <- allocation$outer_partition == "development" & allocation$inner_fold == fold
  training_raw <- candidate[training_index, , drop = FALSE]
  validation_raw <- candidate[validation_index, , drop = FALSE]
  recipe <- r211f_fit_training_recipe(training_raw, "MIMIC", 0.10, 500L)
  training_model <- r211f_apply_fixed_recipe(training_raw, recipe, "MIMIC")
  validation_model <- r211f_apply_fixed_recipe(validation_raw, recipe, "MIMIC")
  predictors <- rfag_predictors(training_model)
  factor_spec <- rfag_factor_spec(recipe)
  training_frame <- rfp_transform_model_frame(rfag_prepare_frame(training_model, predictors, factor_spec))
  validation_frame <- rfp_transform_model_frame(rfag_prepare_frame(validation_model, predictors, factor_spec))
  if (!"paired_mean_saturation_c90__model" %in% predictors ||
      any(grepl("^(spo2_percent|sao2_percent)(__|$)", predictors))) stop("RF_PREDICTOR_CONTRACT", call. = FALSE)
  if (sum(recipe$retained & recipe$feature_role == "candidate") < 1L) stop("RF_REAL_THRESHOLD_RETAINED_NO_CLINICAL_VARIABLE", call. = FALSE)
  for (grid_index in seq_len(nrow(grid))) {
    k <- k + 1L
    mtry <- rff_resolve_mtry(length(predictors), grid$mtry_rule[[grid_index]])
    model <- rff_fit_model(training_frame, mtry, grid$min_node_size[[grid_index]],
                           20260916L + grid_index * 100L + fold, 20L, 1L)
    predicted <- rff_predict_model(model, validation_frame, 1L)
    metric <- rff_metrics(validation_frame$bias_spo2_minus_sao2, predicted)
    fold_rows[[k]] <- data.frame(cohort = "MIMIC", grid_index = grid_index,
      inner_validation_fold = fold, mtry_rule = grid$mtry_rule[[grid_index]],
      min_node_size = grid$min_node_size[[grid_index]], predictor_n = length(predictors),
      resolved_mtry = mtry, rmse = metric[["rmse"]], mae = metric[["mae"]], stringsAsFactors = FALSE)
  }
}
fold_metrics <- do.call(rbind, fold_rows)
if (nrow(fold_metrics) != 60L || any(!is.finite(c(fold_metrics$rmse, fold_metrics$mae)))) stop("RF_INNER_CV_INVALID", call. = FALSE)
tuning <- rff_summarize_tuning(fold_metrics, "MIMIC")
token <- rff_outer_access_token(fold_metrics, tuning$selected, "MIMIC", "synthetic_input", spec$observed_sha256)
expect_error <- function(expr, code) {
  observed <- tryCatch({force(expr); "NO_ERROR"}, error = function(e) conditionMessage(e))
  identical(observed, code)
}
if (!expect_error(rff_prepare_outer(candidate, allocation, data.frame(), "MIMIC", token, "wrong_token"),
                  "RFF_OUTER_ACCESS_GUARD_FAILED")) stop("RF_OUTER_TOKEN_NEGATIVE_TEST_FAILED", call. = FALSE)
allocation_bad <- allocation
allocation_bad$outer_partition[[1L]] <- "test"
if (!expect_error(rff_validate_allocation(candidate, allocation_bad, "MIMIC"),
                  "RFF_OUTER_PATIENT_LEAKAGE")) stop("RF_PATIENT_LEAKAGE_NEGATIVE_TEST_FAILED", call. = FALSE)
if (!expect_error(rfp_round_c90_to_whole_percentage_point(c(1, NA_real_)),
                  "RFP_PAIRED_MEAN_INVALID")) stop("RF_INVALID_OXYGEN_NEGATIVE_TEST_FAILED", call. = FALSE)
outer_raw <- candidate[allocation$outer_partition == "development", , drop = FALSE]
outer_recipe <- r211f_fit_training_recipe(outer_raw, "MIMIC", 0.10, 500L)
recipe_bundle <- cbind(recipe_context = "outer_development", inner_validation_fold = NA_integer_,
                       training_pair_n = nrow(outer_raw), outer_recipe)
outer <- rff_prepare_outer(candidate, allocation, recipe_bundle, "MIMIC", token, token)
outer$development <- rfp_transform_model_frame(outer$development)
outer$test <- rfp_transform_model_frame(outer$test)
outer$test_source_features <- rfp_transform_source_features(outer$test_source_features)
predictors <- outer$predictors
if (sum(outer_recipe$retained & outer_recipe$feature_role == "candidate") < 1L) stop("RF_OUTER_RECIPE_NO_CLINICAL_VARIABLE", call. = FALSE)
if (!"paired_mean_saturation_c90__model" %in% predictors ||
    any(grepl("^(spo2_percent|sao2_percent)(__|$)", predictors))) stop("RF_OUTER_PREDICTOR_CONTRACT", call. = FALSE)
if (!identical(as.numeric(outer$development$paired_mean_saturation_c90__model),
               as.numeric(floor(outer_raw$paired_mean_saturation_c90 + 90 + 0.5) - 90)))
  stop("RF_INTEGER_OXYGEN_EXPECTATION_MISMATCH", call. = FALSE)
if (any(outer$development$paired_mean_saturation_c90__model %% 1 != 0) ||
    any(outer$test$paired_mean_saturation_c90__model %% 1 != 0)) stop("RF_INTEGER_OXYGEN_FAILED", call. = FALSE)
final_mtry <- rff_resolve_mtry(length(predictors), tuning$selected$mtry_rule[[1L]])
final_model <- rff_fit_model(outer$development, final_mtry,
  tuning$selected$min_node_size[[1L]], 20260917L, 50L, 1L)
predicted <- rff_predict_model(final_model, outer$test, 1L)
performance <- rff_metrics(outer$test$bias_spo2_minus_sao2, predicted)
test_meta <- outer$test_meta[, c("analysis_patient_key", "pair_uid")]
bootstrap <- rff_cluster_bootstrap(outer$test$bias_spo2_minus_sao2, predicted,
  test_meta$analysis_patient_key, 50L, 20260918L)
shap_index <- rff_patient_balanced_shap_indices(test_meta$analysis_patient_key,
  test_meta$pair_uid, 15L, 20260919L)
set.seed(20260920L)
shap <- as.matrix(fastshap::explain(final_model,
  X = outer$test[shap_index, predictors, drop = FALSE],
  pred_wrapper = function(object, newdata) rff_predict_model(object, newdata, 1L),
  nsim = 2L, adjust = TRUE, parallel = FALSE))
grouped <- rff_group_shap(shap)
if (length(performance) != 5L || any(!is.finite(performance)) ||
    any(bootstrap$summary$bootstrap_successful != 50L) ||
    length(unique(test_meta$analysis_patient_key[shap_index])) != length(shap_index) ||
    max(abs(rowSums(grouped$values) - rowSums(shap))) > 1e-10) stop("RF_OUTER_QA_FAILURE", call. = FALSE)

checks <- data.frame(check = c("synthetic_patient_split", "inner_5x12_complete",
  "training_only_real_500_pair_retention", "outer_access_token", "no_standalone_oxygen",
  "integer_paired_mean", "outer_test_evaluated_once", "bootstrap_50_complete",
  "patient_balanced_shap", "grouped_shap_additivity", "outer_token_rejection",
  "patient_leakage_rejection", "invalid_oxygen_rejection"), pass = rep(TRUE, 13L))
aggregate <- data.frame(synthetic_patient_n = length(patient), synthetic_pair_n = nrow(candidate),
  development_patient_n = length(development_patient), development_pair_n = nrow(outer$development),
  test_patient_n = length(unique(test_meta$analysis_patient_key)), test_pair_n = nrow(outer$test),
  inner_task_n = nrow(fold_metrics), outer_retained_clinical_n = sum(outer_recipe$retained & outer_recipe$feature_role == "candidate"),
  predictor_n = length(predictors), selected_grid_index = tuning$selected$grid_index[[1L]],
  test_rmse = performance[["rmse"]], test_mae = performance[["mae"]],
  bootstrap_success_n = unique(bootstrap$summary$bootstrap_successful)[[1L]],
  shap_patient_n = length(shap_index), shap_group_n = ncol(grouped$values), stringsAsFactors = FALSE)
write.table(spec, file.path(output_root, "source_receipt.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(checks, file.path(output_root, "checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(aggregate, file.path(output_root, "aggregate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("SYNTHETIC_RF_PASS", file.path(output_root, "status.txt"))
cat("SYNTHETIC_RF_PASS\n")
