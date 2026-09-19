args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: OUTPUT_ROOT", call. = FALSE)
output_root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
if (file.exists(output_root) || dir.exists(output_root)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
extra_library <- Sys.getenv("R9_EXTRA_R_LIB", unset = "")
if (nzchar(extra_library)) {
  if (!dir.exists(extra_library)) stop("EXTRA_R_LIBRARY_MISSING", call. = FALSE)
  .libPaths(c(extra_library, .libPaths()))
}
suppressPackageStartupMessages({
  library(digest)
  library(fastshap)
  library(ranger)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]])
package_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(package_root, "code", "formal", "stage405_rf_amended_formal_core.R"))
source(file.path(package_root, "code", "formal", "stage405_rf_amended_fitting_core.R"))
source(file.path(package_root, "code", "cores", "r211f_core.R"))

set.seed(20260908L)
patient <- rep(sprintf("synthetic_patient_%03d", seq_len(80L)), each = 2L)
n <- length(patient)
stay <- patient
spo2 <- pmin(100, pmax(70, round(stats::rnorm(n, 94, 4), 1)))
sao2 <- pmin(100, pmax(70, spo2 - stats::rnorm(n, 1.5, 2)))
candidate <- data.frame(
  cohort = "MIMIC",
  analysis_patient_key = patient,
  analysis_stay_key = stay,
  pair_uid = sprintf("synthetic_pair_%04d", seq_len(n)),
  bias_spo2_minus_sao2 = spo2 - sao2,
  spo2_percent = spo2,
  sao2_percent = sao2,
  paired_mean_saturation_c90 = ((spo2 + sao2) / 2) - 90,
  age_value_per_10y = rep(seq(3, 9, length.out = 80L), each = 2L),
  age_topcoded_indicator = rep(c(rep(0, 75), rep(1, 5)), each = 2L),
  sex_class = rep(c("female", "male"), length.out = n),
  cci_value = rep(seq(0, 8, length.out = 80L), each = 2L),
  sofa_total = pmax(0, stats::rnorm(n, 8, 3)),
  pao2_candidate_value = pmax(20, stats::rnorm(n, 90, 20)),
  paco2_candidate_value = pmax(15, stats::rnorm(n, 42, 8)),
  ph_candidate_value = pmin(7.8, pmax(6.8, stats::rnorm(n, 7.36, 0.08))),
  lactate_candidate_value = pmax(0.2, stats::rlnorm(n, log(2), 0.4)),
  hb_candidate_value = pmax(4, stats::rnorm(n, 10, 1.5)),
  concurrent_vasoactive_medication_use_state = rep(c("yes", "no"), length.out = n),
  current_invasive_mechanical_ventilation_state = rep(c("no", "yes", "no"), length.out = n),
  stringsAsFactors = FALSE
)
allocation_patient <- unique(patient)
development_patients <- allocation_patient[seq_len(56L)]
allocation <- candidate[, c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid")]
allocation$outer_partition <- ifelse(allocation$analysis_patient_key %in% development_patients, "development", "test")
fold_map <- setNames(rep(seq_len(5L), length.out = length(development_patients)), development_patients)
allocation$inner_fold <- ifelse(allocation$outer_partition == "development", fold_map[allocation$analysis_patient_key], NA_integer_)
development_raw <- candidate[allocation$outer_partition == "development", , drop = FALSE]
recipe <- r211f_fit_training_recipe(development_raw, "MIMIC", 0.10, 500L)

# The formal production threshold is 500 observed pairs. The synthetic fixture is
# deliberately smaller, so only the synthetic recipe's retained flag is changed;
# the production core and production contract are not changed.
recipe$retained[recipe$feature_role == "candidate"] <- TRUE
recipe$observed_n[recipe$retained] <- pmax(recipe$observed_n[recipe$retained], 500L)
for (i in which(recipe$feature_role == "candidate")) {
  feature <- recipe$feature[[i]]
  if (recipe$feature_type[[i]] == "continuous") recipe$fill_value[[i]] <- as.character(stats::median(development_raw[[feature]], na.rm = TRUE))
  if (recipe$feature_type[[i]] == "categorical") {
    value <- r211f_normalize_categorical(development_raw[[feature]])
    recipe$fill_value[[i]] <- r211f_mode(value)
    recipe$allowed_levels[[i]] <- paste(sort(unique(value)), collapse = "|")
  }
}
attr(recipe, "cohort") <- "MIMIC"
attr(recipe, "minimum_support") <- 0.10
attr(recipe, "minimum_observed_n") <- 500L
recipe_bundle <- cbind(recipe_context = "outer_development", inner_validation_fold = NA_integer_, training_pair_n = nrow(development_raw), recipe)

grid <- rff_tuning_grid()
fold_rows <- list()
k <- 0L
for (fold in seq_len(5L)) {
  train_index <- allocation$outer_partition == "development" & allocation$inner_fold != fold
  validation_index <- allocation$outer_partition == "development" & allocation$inner_fold == fold
  training_model <- r211f_apply_fixed_recipe(candidate[train_index, , drop = FALSE], recipe, "MIMIC")
  validation_model <- r211f_apply_fixed_recipe(candidate[validation_index, , drop = FALSE], recipe, "MIMIC")
  predictors <- rfag_predictors(training_model)
  spec <- rfag_factor_spec(recipe)
  training_frame <- rfag_prepare_frame(training_model, predictors, spec)
  validation_frame <- rfag_prepare_frame(validation_model, predictors, spec)
  for (grid_index in seq_len(nrow(grid))) {
    k <- k + 1L
    mtry <- rff_resolve_mtry(length(predictors), grid$mtry_rule[[grid_index]])
    model <- rff_fit_model(training_frame, mtry, grid$min_node_size[[grid_index]], 20260829L + grid_index * 100L + fold, 20L, 1L)
    prediction <- rff_predict_model(model, validation_frame, 1L)
    metric <- rff_metrics(validation_frame$bias_spo2_minus_sao2, prediction)
    fold_rows[[k]] <- data.frame(cohort = "MIMIC", grid_index = grid_index, inner_validation_fold = fold, mtry_rule = grid$mtry_rule[[grid_index]], min_node_size = grid$min_node_size[[grid_index]], predictor_n = length(predictors), resolved_mtry = mtry, rmse = metric[["rmse"]], mae = metric[["mae"]], stringsAsFactors = FALSE)
  }
}
fold_metrics <- do.call(rbind, fold_rows)
tuning <- rff_summarize_tuning(fold_metrics, "MIMIC")
token <- rff_outer_access_token(fold_metrics, tuning$selected, "MIMIC", "synthetic_input", "synthetic_code")
outer <- rff_prepare_outer(candidate, allocation, recipe_bundle, "MIMIC", token, token)
predictors <- outer$predictors
development_frame <- outer$development
test_frame <- outer$test
final_mtry <- rff_resolve_mtry(length(predictors), tuning$selected$mtry_rule[[1L]])
final_model <- rff_fit_model(development_frame, final_mtry, tuning$selected$min_node_size[[1L]], 20260932L, 50L, 1L)
prediction <- rff_predict_model(final_model, test_frame, 1L)
test_meta <- outer$test_meta[, c("analysis_patient_key", "pair_uid")]
bootstrap <- rff_cluster_bootstrap(test_frame$bias_spo2_minus_sao2, prediction, test_meta$analysis_patient_key, 50L, 20261032L)
shap_index <- rff_patient_balanced_shap_indices(test_meta$analysis_patient_key, test_meta$pair_uid, 1000L, 20261132L)
set.seed(20261232L)
shap <- as.matrix(fastshap::explain(final_model, X = test_frame[shap_index, predictors, drop = FALSE], pred_wrapper = function(object, newdata) rff_predict_model(object, newdata, 1L), nsim = 3L, adjust = TRUE, parallel = FALSE))
grouped <- rff_group_shap(shap)

checks <- data.frame(
  check = c("tuning_60_of_60", "outer_token_created", "no_standalone_spo2", "no_standalone_sao2", "paired_mean_present", "performance_five_metrics", "bootstrap_50_complete", "shap_one_pair_per_patient", "grouped_shap_additive"),
  status = c(
    if (nrow(fold_metrics) == 60L) "PASS" else "FAIL",
    if (nchar(token) == 64L) "PASS" else "FAIL",
    if (!any(grepl("^spo2_percent", predictors))) "PASS" else "FAIL",
    if (!any(grepl("^sao2_percent", predictors))) "PASS" else "FAIL",
    if ("paired_mean_saturation_c90__model" %in% predictors) "PASS" else "FAIL",
    if (length(rff_metrics(test_frame$bias_spo2_minus_sao2, prediction)) == 5L) "PASS" else "FAIL",
    if (all(bootstrap$summary$bootstrap_successful == 50L)) "PASS" else "FAIL",
    if (length(unique(test_meta$analysis_patient_key[shap_index])) == length(shap_index)) "PASS" else "FAIL",
    if (max(abs(rowSums(grouped$values) - rowSums(shap))) < 1e-10) "PASS" else "FAIL"
  ),
  stringsAsFactors = FALSE
)
utils::write.table(checks, file.path(output_root, "end_to_end_synthetic_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (any(checks$status != "PASS")) stop("RFF_END_TO_END_SYNTHETIC_HOLD", call. = FALSE)
cat("RFF_END_TO_END_SYNTHETIC_PASS\n")
