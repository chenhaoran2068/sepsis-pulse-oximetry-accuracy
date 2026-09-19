args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: COHORT WINDOW_MINUTES ANALYSIS_INPUT_DIR NEW_OUTPUT_DIR", call. = FALSE)
cohort <- args[[1L]]
window <- suppressWarnings(as.integer(args[[2L]]))
input <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[4L]], winslash = "/", mustWork = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (!cohort %in% cohorts || is.na(window) || !window %in% c(60L, 5L))
  stop("FORMAL_MODEL_COHORT_OR_WINDOW_INVALID", call. = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (!identical(readLines(file.path(input, "status.txt"), warn = FALSE),
               "ANALYSIS_INPUT_BUILD_PASS")) stop("ANALYSIS_INPUT_GATE_NOT_PASS", call. = FALSE)
input_path <- file.path(input, paste0("lmm_input_", window, ".rds"))
if (!file.exists(input_path)) stop("LMM_INPUT_MISSING", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)
extra_library <- Sys.getenv("R9_EXTRA_R_LIB", unset = "")
if (nzchar(extra_library)) {
  if (!dir.exists(extra_library)) stop("EXTRA_R_LIBRARY_MISSING", call. = FALSE)
  .libPaths(c(extra_library, .libPaths()))
}
suppressPackageStartupMessages({
  library(data.table); library(digest); library(mice); library(miceadds)
  library(lme4); library(sandwich); library(car); library(posterior)
})
for (name in c("stage405_multilevel_core.R", "stage405_formal_convergence_adapter.R",
               "stage405_formal_imputation.R", "stage405_formal_lmm_stream.R"))
  source(file.path(root, "code", "formal", name), local = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (name in c("aggregate", "imputations", "checkpoints", "runtime"))
  dir.create(file.path(output, name), showWarnings = FALSE)
writeLines("FORMAL_MI_LMM_CANDIDATE_RUNNING", file.path(output, "status.txt"))
write_table <- function(x, name) data.table::fwrite(x, file.path(output, "aggregate", name), sep = "\t")
sha <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
tryCatch({
  x <- data.table::as.data.table(readRDS(input_path))
  if (any(x$cohort != cohort) || any(x$pair_window_minutes != window) || !nrow(x))
    stop("FORMAL_MODEL_INPUT_SCOPE_INVALID", call. = FALSE)
  coverage <- s405m_coverage_ledger(x, cohort)
  predictors <- s405m_select_predictors(coverage)
  if (!length(predictors)) stop("FORMAL_MODEL_NO_ELIGIBLE_PREDICTOR", call. = FALSE)
  seed <- 202609070L + match(cohort, cohorts) * 100L + match(window, c(60L, 5L))
  warnings <- character()
  imputed <- withCallingHandlers(
    s405f_run_imputation(x, cohort, predictors, 30L, 20L, seed, FALSE),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  validation <- s405m_validate_imputation(x, imputed, cohort)
  domain <- s405m_imputed_domain_audit(imputed)
  if (length(imputed$completed) != 30L ||
      !all(validation$variable$observed_unchanged) ||
      !all(validation$variable$completed_no_missing) ||
      !all(validation$stay_constancy$stay_constant) ||
      (nrow(validation$patient_sex_constancy) && !all(validation$patient_sex_constancy$patient_constant)) ||
      !all(domain$domain_pass) ||
      any(imputed$convergence$gate$blocking_failure_n > 0L))
    stop("FORMAL_IMPUTATION_QA_BLOCK", call. = FALSE)
  write_table(coverage, "variable_coverage.tsv")
  write_table(validation$variable, "imputation_variable_validation.tsv")
  write_table(validation$stay_constancy, "stay_value_constancy.tsv")
  write_table(validation$patient_sex_constancy, "patient_sex_constancy.tsv")
  write_table(domain, "imputed_value_domain_audit.tsv")
  write_table(imputed$convergence$summary, "convergence_summary.tsv")
  write_table(imputed$convergence$gate, "convergence_gate.tsv")
  write_table(data.table(message = unique(warnings),
                         count = tabulate(match(warnings, unique(warnings)),
                                          nbins = length(unique(warnings)))),
              "runtime_warnings.tsv")
  output_rows <- vector("list", 30L)
  fit_results <- vector("list", 30L)
  for (i in seq_len(30L)) {
    completed <- imputed$completed[[i]]
    dataset_path <- file.path(output, "imputations", sprintf("completed_%02d.rds", i))
    saveRDS(completed, dataset_path, version = 3)
    output_rows[[i]] <- data.frame(imputation = i, file = basename(dataset_path),
                                   row_n = nrow(completed), sha256 = sha(dataset_path))
    fit_results[[i]] <- list(
      primary = s405l_fit_primary_one(completed, cohort, predictors),
      balanced = s405l_fit_balanced_one(completed, predictors))
    saveRDS(fit_results[[i]], file.path(output, "checkpoints", sprintf("fit_%02d.rds", i)),
            version = 3)
  }
  write_table(data.table::rbindlist(output_rows), "completed_dataset_manifest.tsv")
  patient_n <- data.table::uniqueN(imputed$completed[[1L]]$cluster_patient)
  pooled <- s405l_pool_stream_results(fit_results, patient_n)
  collect <- function(kind, part) data.table::rbindlist(lapply(seq_along(fit_results), function(i) {
    x <- data.table::as.data.table(fit_results[[i]][[kind]][[part]])
    x[, imputation := i]
    x
  }), fill = TRUE)
  fit_status <- collect("primary", "status")
  design <- collect("primary", "design")
  gvif <- collect("primary", "gvif")
  residual <- collect("primary", "residual")
  weight_status <- collect("balanced", "weight_status")
  fit_checks <- data.table::data.table(
    check = c("exact_30_fits", "all_lmms_converged", "all_hessians_positive_definite",
              "no_singular_fit", "random_variance_floor_pass", "all_designs_full_rank",
              "all_adjusted_gvif_below_5", "finite_residuals_and_fitted",
              "balanced_patient_total_weights_one", "primary_mcse_ratio_at_most_0_10",
              "balanced_mcse_ratio_at_most_0_10"),
    pass = c(nrow(fit_status) == 30L, all(fit_status$converged),
             all(fit_status$hessian_positive_definite), !any(fit_status$singular),
             all(fit_status$random_variance_floor_pass), all(design$full_rank),
             all(gvif$adjusted_gvif < 5),
             all(residual$all_residuals_finite & residual$all_fitted_finite),
             all(weight_status$all_patient_total_weights_equal_one),
             all(pooled$primary$monte_carlo_error_to_se <= 0.10),
             all(pooled$balanced$monte_carlo_error_to_se <= 0.10)))
  write_table(fit_status, "primary_fit_status.tsv")
  write_table(design, "primary_design_status.tsv")
  write_table(gvif, "primary_gvif.tsv")
  write_table(residual, "primary_residual_status.tsv")
  write_table(weight_status, "patient_balanced_weight_status.tsv")
  write_table(fit_checks, "formal_lmm_gate.tsv")
  if (!all(fit_checks$pass))
    stop(paste("FORMAL_LMM_GATE_FAILED:", paste(fit_checks$check[!fit_checks$pass], collapse = ",")), call. = FALSE)
  write_table(s405l_apply_reporting_scale(pooled$primary), "lmm_primary_pooled.tsv")
  write_table(s405l_apply_reporting_scale(pooled$balanced), "lmm_patient_balanced_pooled.tsv")
  write_table(data.table(cohort = cohort, window_minutes = window, pair_n = nrow(x),
                         patient_n = patient_n, imputations = 30L, iterations = 20L,
                         seed = seed, input_sha256 = sha(input_path),
                         convergence_review_required_n =
                           imputed$convergence$gate$review_required_n[[1L]],
                         candidate_only = TRUE), "run_summary.tsv")
  writeLines("FORMAL_MI_LMM_CANDIDATE_COMPLETE_QA_PENDING",
             file.path(output, "status.txt"))
  cat("FORMAL_MI_LMM_CANDIDATE_COMPLETE_QA_PENDING\n")
}, error = function(e) {
  writeLines(paste("FORMAL_MI_LMM_CANDIDATE_FAIL", conditionMessage(e)),
             file.path(output, "status.txt"))
  stop(e)
})
