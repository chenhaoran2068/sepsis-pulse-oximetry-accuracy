args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L || !args[[5L]] %in% c("paper", "smoke"))
  stop("Usage: COHORT WINDOW_MINUTES ANALYSIS_INPUT_DIR NEW_OUTPUT_DIR paper|smoke", call. = FALSE)
cohort <- args[[1L]]
window <- as.integer(args[[2L]])
input <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[4L]], winslash = "/", mustWork = FALSE)
mode <- args[[5L]]
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (!cohort %in% cohorts || !window %in% c(60L, 5L)) stop("AGREEMENT_SCOPE_INVALID", call. = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)
source(file.path(root, "code", "cores", "agreement_core.R"))
source(file.path(root, "code", "cores", "stage404_core.R"))
suppressPackageStartupMessages({library(data.table); library(digest); library(lme4); library(nlme)})
pair_path <- file.path(input, paste0("matched_pairs_", window, ".rds"))
analytic_path <- file.path(input, paste0("analytic_pairs_", window, ".rds"))
receipt <- read.delim(file.path(input, "receipt.tsv"), check.names = FALSE)
hash_ok <- function(path) {
  row <- receipt[receipt$file == basename(path), ]
  file.exists(path) && nrow(row) == 1L &&
    tolower(row$sha256[[1L]]) == tolower(digest::digest(file = path, algo = "sha256"))
}
if (!hash_ok(pair_path) || !hash_ok(analytic_path))
  stop("AGREEMENT_PAIR_HASH_INVALID", call. = FALSE)
x <- data.table::as.data.table(readRDS(pair_path))
analytic <- data.table::as.data.table(readRDS(analytic_path))
if (nrow(x) != nrow(analytic) ||
    !identical(as.character(x$stay_key_internal), as.character(analytic$stay_key_internal)) ||
    !identical(as.character(x$patient_key_internal), as.character(analytic$patient_key_internal)))
  stop("AGREEMENT_ANALYTIC_PAIR_IDENTITY_MISMATCH", call. = FALSE)
x[, pair_key_internal := analytic$pair_key_internal]
validate_pair_input(x, cohort, if (window == 60L) "M60" else "M5")
z <- prepare_agreement_data(x)
nested <- any(z[, .(stay_n = uniqueN(stay)), by = patient]$stay_n > 1L)
if (nested != (cohort == "eICU")) stop("AGREEMENT_RANDOM_STRUCTURE_MISMATCH", call. = FALSE)
reps <- if (mode == "paper") 2000L else 30L
cohort_index <- match(cohort, cohorts)
base_seed <- 403000L + (if (window == 60L) 0L else 10000L) + cohort_index * 100L
dir.create(output, recursive = TRUE)
writeLines("AGREEMENT_RUNNING", file.path(output, "status.txt"))
write_tsv <- function(x, name) data.table::fwrite(x, file.path(output, name), sep = "\t")
tryCatch({
  overall <- summarize_overall(z, b = reps, seed = base_seed + 1L, nested = nested)
  overall_row <- data.table(cohort = cohort, window_minutes = window, mode = mode,
    pair_n = nrow(z), patient_n = uniqueN(z$patient), stay_n = uniqueN(z$stay),
    mean_bias = overall$fit$mean_bias,
    mean_bias_ci_lower = overall$ci$mean_bias[[1L]],
    mean_bias_ci_upper = overall$ci$mean_bias[[2L]],
    between_patient_variance = overall$fit$between_patient_variance,
    between_stay_within_patient_variance = overall$fit$between_stay_within_patient_variance,
    within_stay_residual_variance = overall$fit$within_stay_residual_variance,
    total_difference_sd = overall$fit$total_sd,
    lower_loa = overall$fit$lower_loa,
    lower_loa_ci_lower = overall$ci$lower_loa[[1L]],
    lower_loa_ci_upper = overall$ci$lower_loa[[2L]],
    upper_loa = overall$fit$upper_loa,
    upper_loa_ci_lower = overall$ci$upper_loa[[1L]],
    upper_loa_ci_upper = overall$ci$upper_loa[[2L]],
    arms = overall$dist$pair_arms,
    arms_ci_lower = overall$ci$pair_arms[[1L]],
    arms_ci_upper = overall$ci$pair_arms[[2L]],
    empirical_p2_5 = overall$dist$empirical_p2_5,
    empirical_p2_5_ci_lower = overall$ci$empirical_p2_5[[1L]],
    empirical_p2_5_ci_upper = overall$ci$empirical_p2_5[[2L]],
    empirical_p97_5 = overall$dist$empirical_p97_5,
    empirical_p97_5_ci_lower = overall$ci$empirical_p97_5[[1L]],
    empirical_p97_5_ci_upper = overall$ci$empirical_p97_5[[2L]],
    bootstrap_requested = reps, bootstrap_successful = sum(overall$boot$success))
  write_tsv(overall_row, "agreement_overall.tsv")
  write_tsv(data.table(
    cohort = cohort, window_minutes = window, mode = mode,
    patient_n = uniqueN(z$patient), pair_n = nrow(z),
    balanced_mean_bias = overall$dist$balanced_mean,
    balanced_mean_bias_ci_lower = overall$ci$balanced_mean[[1L]],
    balanced_mean_bias_ci_upper = overall$ci$balanced_mean[[2L]],
    balanced_lower_loa = overall$dist$balanced_lower_loa,
    balanced_lower_loa_ci_lower = overall$ci$balanced_lower_loa[[1L]],
    balanced_lower_loa_ci_upper = overall$ci$balanced_lower_loa[[2L]],
    balanced_upper_loa = overall$dist$balanced_upper_loa,
    balanced_upper_loa_ci_lower = overall$ci$balanced_upper_loa[[1L]],
    balanced_upper_loa_ci_upper = overall$ci$balanced_upper_loa[[2L]],
    balanced_arms = overall$dist$balanced_arms,
    balanced_arms_ci_lower = overall$ci$balanced_arms[[1L]],
    balanced_arms_ci_upper = overall$ci$balanced_arms[[2L]],
    balanced_empirical_p2_5 = overall$dist$balanced_empirical_p2_5,
    balanced_empirical_p97_5 = overall$dist$balanced_empirical_p97_5,
    bootstrap_requested = reps, bootstrap_successful = sum(overall$boot$success)
  ), "agreement_patient_balanced.tsv")
  strata <- lapply(seq_along(s404_strata()), function(i) {
    name <- s404_strata()[[i]]
    subset <- z[spo2_stratum == name]
    if (nrow(subset) == 0L || uniqueN(subset$patient) < 2L)
      stop(paste("AGREEMENT_STRATUM_NOT_ESTIMABLE", name), call. = FALSE)
    result <- summarize_stratum(subset, b = reps, seed = base_seed + 10L + i, nested = nested)
    data.table(cohort = cohort, window_minutes = window, spo2_stratum = name,
      pair_n = nrow(subset), patient_n = uniqueN(subset$patient),
      mean_bias = result$fit$mean_bias,
      mean_bias_ci_lower = result$ci[[1L]][[1L]],
      mean_bias_ci_upper = result$ci[[1L]][[2L]],
      lower_loa = result$fit$lower_loa,
      lower_loa_ci_lower = result$ci[[2L]][[1L]],
      lower_loa_ci_upper = result$ci[[2L]][[2L]],
      upper_loa = result$fit$upper_loa,
      upper_loa_ci_lower = result$ci[[3L]][[1L]],
      upper_loa_ci_upper = result$ci[[3L]][[2L]],
      bootstrap_successful = sum(result$boot$success))
  })
  write_tsv(rbindlist(strata), "agreement_by_spo2_stratum.tsv")
  proportional <- fit_proportional_bias(z, nested = nested)
  write_tsv(data.table(cohort = cohort, window_minutes = window,
                       slope = proportional$slope, ci_lower = proportional$lower,
                       ci_upper = proportional$upper,
                       model_estimated_bias_at_90 = proportional$intercept_at_90,
                       singular_fit = proportional$singular_fit), "agreement_proportional_bias.tsv")
  hetero <- fit_heteroscedasticity(z, nested = nested)
  write_tsv(data.table(cohort = cohort, window_minutes = window,
                       variance_function_parameter = hetero$delta,
                       ci_lower = hetero$delta_lower, ci_upper = hetero$delta_upper,
                       residual_sd_at_90 = hetero$residual_sigma_at_90,
                       likelihood_ratio = hetero$likelihood_ratio, p_value = hetero$p_value),
            "agreement_heteroscedasticity.tsv")
  write_tsv(data.table(cohort = cohort, window_minutes = window, mode = mode,
                       input_sha256 = digest::digest(file = pair_path, algo = "sha256"),
                       nested_random_effects = nested,
                       bootstrap_replicates = reps), "run_receipt.tsv")
  writeLines(if (mode == "paper") "AGREEMENT_PAPER_SETTINGS_COMPLETE_CANDIDATE_ONLY" else
               "AGREEMENT_SMOKE_COMPLETE_NOT_PAPER_SETTINGS", file.path(output, "status.txt"))
  cat("AGREEMENT_RUN_PASS mode=", mode, "\n")
}, error = function(e) {
  writeLines(paste("AGREEMENT_FAIL", conditionMessage(e)), file.path(output, "status.txt"))
  stop(e)
})
