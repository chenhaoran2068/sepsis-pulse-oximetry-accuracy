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
writeLines("SYNTHETIC_MI_LMM_RUNNING", file.path(output_root, "status.txt"))
options(error = function() {
  writeLines("SYNTHETIC_MI_LMM_FAIL", file.path(output_root, "status.txt"))
  quit(save = "no", status = 1L)
})

suppressPackageStartupMessages({
  library(data.table); library(mice); library(miceadds); library(lme4); library(sandwich); library(car)
})
module_root <- file.path(candidate_root, "code", "cores")
spec <- data.frame(module = c("R211B", "MI_PILOT", "LMM"), relative_path = c(
  "r211b_core.R",
  "stage405_mi_demo_core.R",
  "stage405_lmm_demo_core.R"
), expected_sha256 = c(
  "68F51663252389676FBC58BAF7CE34B6671F155D9A7801C66B47A829EEE99C09",
  "0D02CC6810CA8CDECF67F6F91575CE4D1FE342B5B05803C18B19827005D81167",
  "A2AB3BFA453378ECEB0D373A21E84D32B3CBE34FFA70A8090341802B35E21BDC"
), stringsAsFactors = FALSE)
paths <- file.path(module_root, spec$relative_path)
if (any(!file.exists(paths))) stop("MODULE_SOURCE_MISSING", call. = FALSE)
spec$observed_sha256 <- toupper(unname(tools::sha256sum(paths)))
if (any(spec$observed_sha256 != spec$expected_sha256)) stop("MODULE_SOURCE_HASH_MISMATCH", call. = FALSE)
for (path in paths) source(path, local = FALSE)
source(file.path(candidate_root, "src", "synthetic_model_inputs.R"), local = FALSE)
original_two_level_imputer <- mice.impute.2l.pmm.scaled
singular_by_missing_n <- integer()
mice.impute.2l.pmm.scaled <- function(y, ry, x, type, ...) {
  missing_n <- as.character(sum(!ry))
  withCallingHandlers(original_two_level_imputer(y = y, ry = ry, x = x, type = type, ...),
    message = function(msg) {
      if (grepl("boundary \\(singular\\) fit", conditionMessage(msg))) {
        if (!missing_n %in% names(singular_by_missing_n)) singular_by_missing_n[[missing_n]] <<- 0L
        singular_by_missing_n[[missing_n]] <<- singular_by_missing_n[[missing_n]] + 1L
      }
    })
}

upstream <- readRDS(file.path(upstream_root, "synthetic_pairs_INTERNAL_ONLY.rds"))
if (!is.list(upstream) || !all(c("M60", "M5") %in% names(upstream))) stop("SYNTHETIC_PAIR_INTERFACE_INVALID", call. = FALSE)
bundle <- private_build_model_inputs(upstream$M60)
all_pairs <- bundle$pairs
all_candidate <- bundle$candidate
if (nrow(all_pairs) != nrow(all_candidate)) stop("R211B_CARDINALITY_CHANGED", call. = FALSE)

patient_ids <- unique(all_candidate$analysis_patient_key)
if (length(patient_ids) < 60L) stop("INSUFFICIENT_SYNTHETIC_PATIENTS", call. = FALSE)
keep <- all_candidate$analysis_patient_key %in% patient_ids[seq_len(60L)]
x <- as.data.table(all_candidate[keep, , drop = FALSE])
source_pair <- all_pairs[match(x$pair_uid, all_pairs$pair_key_internal), , drop = FALSE]
if (anyNA(source_pair$pair_key_internal)) stop("PAIR_TO_MODEL_INPUT_UNMATCHED", call. = FALSE)
x[, `:=`(
  patient_key_internal = analysis_patient_key,
  stay_key_internal = analysis_stay_key,
  sao2_anchor_icu_hours = source_pair$sao2_relative_icu_minutes / 60
)]
x <- s405m_derive_age_fields(x, "MIMIC")
missing_patients <- unique(x$patient_key_internal)
x[patient_key_internal %in% missing_patients[c(3L, 11L, 19L)], age_value_per_10y := NA_real_]
x[patient_key_internal %in% missing_patients[c(6L, 17L)], sex_class := NA_character_]
x[patient_key_internal %in% missing_patients[c(8L, 22L)], cci_value := NA_real_]
x[seq(5L, .N, by = 23L), sofa_total := NA_real_]
x[seq(7L, .N, by = 19L), pao2_candidate_value := NA_real_]
x[seq(9L, .N, by = 17L), paco2_candidate_value := NA_real_]
x[seq(11L, .N, by = 29L), ph_candidate_value := NA_real_]
eligible <- c("age_value_per_10y", "age_topcoded_indicator", "sex_class", "cci_value",
  "sofa_total", "pao2_candidate_value", "paco2_candidate_value", "ph_candidate_value")
imputation_internal_singular_message_n <- 0L
result <- withCallingHandlers(
  s405m_run_imputation(x, "MIMIC", eligible, m = 3L, maxit = 5L,
    seed = 20260916L, print_flag = FALSE, execution_mode = "pilot"),
  message = function(msg) {
    if (grepl("boundary \\(singular\\) fit", conditionMessage(msg), fixed = FALSE))
      imputation_internal_singular_message_n <<- imputation_internal_singular_message_n + 1L
  })
if (nrow(result$logged_events)) {
  cat("SYNTHETIC_MI_LOGGED_EVENT_COUNT=", nrow(result$logged_events), "\n", sep = "")
}
cat("SYNTHETIC_SINGULAR_BY_MISSING_N=", paste(names(singular_by_missing_n), singular_by_missing_n, sep = ":", collapse = ","), "\n", sep = "")
if (nrow(result$logged_events) > 0L) stop("SYNTHETIC_MI_LOGGED_EVENTS", call. = FALSE)
if (imputation_internal_singular_message_n > 0L || sum(singular_by_missing_n) > 0L)
  stop("SYNTHETIC_MI_INTERMEDIATE_SINGULAR", call. = FALSE)
validation <- s405m_validate_imputation(x, result, "MIMIC")
domain <- s405m_imputed_domain_audit(result)
if (length(result$completed) != 3L || result$m != 3L || result$maxit != 5L) stop("PILOT_DESIGN_CHANGED", call. = FALSE)
if (!all(validation$variable$observed_unchanged) || !all(validation$variable$completed_no_missing) ||
    !all(validation$stay_constancy$stay_constant) || !all(domain$domain_pass)) stop("IMPUTATION_QA_FAILURE", call. = FALSE)
if (any(vapply(result$completed, function(d) anyNA(d[, eligible, drop = FALSE]), logical(1L)))) stop("COMPLETED_DATA_MISSING", call. = FALSE)
if (!identical(result$spec$convergence_interpretation, "pilot_interface_screen_only_not_formal_convergence_evidence")) stop("PILOT_MISLABELED_FORMAL", call. = FALSE)

fits <- lapply(result$completed, function(d) list(
  primary = s405l_fit_primary_one(as.data.table(d), "MIMIC", eligible),
  balanced = s405l_fit_balanced_one(as.data.table(d), eligible)
))
if (any(!vapply(fits, function(z) z$primary$status$converged[[1L]], logical(1L)))) stop("SYNTHETIC_LMM_NONCONVERGENCE", call. = FALSE)
singular_n <- sum(vapply(fits, function(z) z$primary$status$singular[[1L]], logical(1L)))
cat("SYNTHETIC_LMM_SINGULAR_COUNT=", singular_n, "\n", sep = "")
if (singular_n > 0L) stop("SYNTHETIC_LMM_SINGULAR", call. = FALSE)
if (any(!vapply(fits, function(z) z$balanced$weight_status$all_patient_total_weights_equal_one[[1L]], logical(1L)))) stop("PATIENT_BALANCE_WEIGHT_FAILURE", call. = FALSE)
pooled <- s405l_pool_stream_results(fits, length(unique(x$patient_key_internal)))
if (any(!is.finite(pooled$primary$estimate)) || any(!is.finite(pooled$balanced$estimate))) stop("POOLED_ESTIMATE_NONFINITE", call. = FALSE)
if (grepl("spo2_saturation_percent|sao2_saturation_percent", paste(deparse(fits[[1L]]$primary$formula), collapse = " "))) stop("STANDALONE_OXYGEN_IN_LMM", call. = FALSE)

checks <- data.frame(check = c("R211B_candidate_cardinality", "pilot_3x5_only", "observed_values_unchanged",
  "all_completed_nonmissing", "stay_level_constancy", "imputed_domains_valid", "LMM_all_three_converged",
  "balanced_patient_weights", "pooled_estimates_finite", "no_standalone_oxygen_in_LMM"),
  pass = c(nrow(all_pairs) == nrow(all_candidate), result$m == 3L && result$maxit == 5L,
    all(validation$variable$observed_unchanged), all(validation$variable$completed_no_missing),
    all(validation$stay_constancy$stay_constant), all(domain$domain_pass),
    all(vapply(fits, function(z) z$primary$status$converged[[1L]], logical(1L))) && singular_n == 0L,
    all(vapply(fits, function(z) z$balanced$weight_status$all_patient_total_weights_equal_one[[1L]], logical(1L))),
    all(is.finite(pooled$primary$estimate)) && all(is.finite(pooled$balanced$estimate)),
    !grepl("spo2_saturation_percent|sao2_saturation_percent", paste(deparse(fits[[1L]]$primary$formula), collapse = " "))),
  stringsAsFactors = FALSE)
aggregate <- data.frame(synthetic_patient_n = length(unique(x$patient_key_internal)),
  synthetic_pair_n = nrow(x), pilot_imputations = result$m, pilot_iterations = result$maxit,
  imputed_cell_n = sum(domain$imputed_cell_n), imputation_domain_fail_n = sum(!domain$domain_pass),
  imputation_logged_event_n = nrow(result$logged_events),
  imputation_internal_singular_message_n = imputation_internal_singular_message_n,
  lmm_singular_fit_n = singular_n,
  lmm_fit_n = length(fits), lmm_term_n = nrow(pooled$primary),
  primary_all_finite = all(is.finite(pooled$primary$estimate)),
  balanced_all_finite = all(is.finite(pooled$balanced$estimate)), stringsAsFactors = FALSE)
write.table(spec, file.path(output_root, "source_receipt.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(checks, file.path(output_root, "checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(aggregate, file.path(output_root, "aggregate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
event_summary <- if (nrow(result$logged_events)) {
  as.data.frame(table(dependent_variable = as.character(result$logged_events$dep),
                      event_method = as.character(result$logged_events$meth),
                      event_detail = as.character(result$logged_events$out)), stringsAsFactors = FALSE)
} else data.frame(dependent_variable = character(), event_method = character(),
                  event_detail = character(), Freq = integer())
event_summary <- event_summary[event_summary$Freq > 0L, , drop = FALSE]
names(event_summary)[names(event_summary) == "Freq"] <- "event_n"
write.table(event_summary, file.path(output_root, "imputation_logged_event_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("SYNTHETIC_MI_LMM_PASS", file.path(output_root, "status.txt"))
cat("SYNTHETIC_MI_LMM_PASS\n")
