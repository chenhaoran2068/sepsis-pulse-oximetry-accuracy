args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: OUTPUT_ROOT", call. = FALSE)
output_root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
dir.create(file.path(output_root, "aggregate"), recursive = TRUE, showWarnings = FALSE)
package_root <- normalizePath(file.path(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][[1L]])), ".."), winslash = "/", mustWork = TRUE)
extra_library <- Sys.getenv("R9_EXTRA_R_LIB", unset = "")
if (nzchar(extra_library)) {
  if (!dir.exists(extra_library)) stop("EXTRA_R_LIBRARY_MISSING", call. = FALSE)
  .libPaths(c(extra_library, .libPaths()))
}
suppressPackageStartupMessages({library(data.table); library(lme4); library(sandwich); library(digest); library(arrow); library(miceadds); library(Matrix); library(car); library(mice)})
core_path <- file.path(package_root, "code", "formal", "stage405_multilevel_core.R")
source(core_path)
source(file.path(package_root, "code", "formal", "stage405_formal_lmm_stream.R"))

set.seed(20260908)
make_data <- function(eicu = FALSE) {
  patient <- rep(seq_len(80), each = if (eicu) 8 else 6)
  stay <- if (eicu) rep(rep(seq_len(160), each = 4), length.out = length(patient)) else patient
  n <- length(patient)
  sex <- factor(ifelse(patient %% 2 == 0, "female", "male"), levels = c("female", "male"))
  pm <- rnorm(n, 3, 2)
  sofa <- pmax(0, round(rnorm(n, 7, 3)))
  pao2 <- pmax(30, rnorm(n, 90, 20))
  ph <- pmin(7.8, pmax(6.8, rnorm(n, 7.36, .08)))
  pe <- rnorm(max(patient), 0, .7)[patient]
  se <- if (eicu) rnorm(max(stay), 0, .4)[stay] else 0
  data.table(
    bias_spo2_minus_sao2 = .3 * pm + .12 * sofa + .2 * (sex == "female") - .008 * pao2 + 1.2 * (ph - 7.4) + pe + se + rnorm(n),
    paired_mean_saturation_c90 = pm, sao2_anchor_icu_hours = runif(n, 0, 168),
    age_value_per_10y = 6 + patient / 100, sex_class = sex, cci_value = patient %% 5,
    sofa_total = sofa, pao2_candidate_value = pao2, ph_candidate_value = ph,
    cluster_stay = stay, cluster_patient = patient
  )
}

checks <- list()
add <- function(name, pass, detail = "") checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail))
for (cohort in c("MIMIC", "eICU")) {
  d <- make_data(cohort == "eICU")
  completed <- lapply(1:3, function(i) copy(d)[, sofa_total := pmax(0, sofa_total + sample(c(-1, 0, 1), .N, TRUE, c(.05, .9, .05)))])
  predictors <- c("age_value_per_10y", "sex_class", "cci_value", "sofa_total", "pao2_candidate_value", "ph_candidate_value")
  streamed <- lapply(completed, function(x) list(primary = s405l_fit_primary_one(x, cohort, predictors), balanced = s405l_fit_balanced_one(x, predictors)))
  pooled <- s405l_pool_stream_results(streamed, uniqueN(d$cluster_patient))
  expected_components <- if (cohort == "eICU") 2L else 1L
  add(paste0(cohort, "_formula_excludes_individual_saturations"), !grepl("spo2_saturation_percent|sao2_saturation_percent", paste(deparse(streamed[[1L]]$primary$formula), collapse = " ")))
  add(paste0(cohort, "_female_vs_male_term"), "sex_classfemale" %in% streamed[[1L]]$primary$terms && !"sex_classmale" %in% streamed[[1L]]$primary$terms)
  add(paste0(cohort, "_random_component_count"), all(vapply(streamed, function(x) x$primary$status$random_component_n[[1L]] == expected_components, logical(1L))))
  add(paste0(cohort, "_fixed_design_full_rank"), all(vapply(streamed, function(x) x$primary$design$full_rank[[1L]], logical(1L))))
  add(paste0(cohort, "_gvif_finite"), all(vapply(streamed, function(x) all(is.finite(x$primary$gvif$adjusted_gvif)), logical(1L))))
  add(paste0(cohort, "_residual_diagnostics_finite"), all(vapply(streamed, function(x) x$primary$residual$all_residuals_finite[[1L]] && x$primary$residual$all_fitted_finite[[1L]], logical(1L))))
  add(paste0(cohort, "_patient_weights_equal_one"), all(vapply(streamed, function(x) x$balanced$weight_status$all_patient_total_weights_equal_one[[1L]], logical(1L))))
  add(paste0(cohort, "_pooled_primary_finite"), all(is.finite(pooled$primary$estimate) & is.finite(pooled$primary$std_error)))
  add(paste0(cohort, "_pooled_balanced_finite"), all(is.finite(pooled$balanced$estimate) & is.finite(pooled$balanced$std_error)))
  test_term <- pooled$primary$term[[2L]]
  x <- rbindlist(lapply(seq_along(streamed), function(i) data.table(imputation = i, estimate = streamed[[i]]$primary$coefficients[[test_term]], variance = streamed[[i]]$primary$variances[[test_term]])))
  independent <- mice::pool.scalar(x$estimate, x$variance, n = uniqueN(d$cluster_patient), k = nrow(pooled$primary))
  candidate_row <- pooled$primary[term == test_term]
  add(paste0(cohort, "_rubin_pool_matches_mice"), isTRUE(all.equal(candidate_row$estimate, independent$qbar, tolerance = 1e-10)) && isTRUE(all.equal(candidate_row$std_error^2, independent$t, tolerance = 1e-10)))
}

scale_input <- data.table(
  term = c("pao2_candidate_value", "paco2_candidate_value", "ph_candidate_value", "age_value_per_10y", "sex_classfemale"),
  estimate = 1, std_error = 2, ci_lower = -3, ci_upper = 4,
  monte_carlo_error_estimate = .5
)
scaled <- s405l_apply_reporting_scale(scale_input)
add("reporting_scale_pao2_per_10", scaled[term == "pao2_candidate_value", estimate] == 10)
add("reporting_scale_paco2_per_10", scaled[term == "paco2_candidate_value", estimate] == 10)
add("reporting_scale_ph_per_0_10", scaled[term == "ph_candidate_value", estimate] == .1)
add("reporting_scale_age_already_per_10", scaled[term == "age_value_per_10y", estimate] == 1)
add("reporting_scale_female_vs_male", scaled[term == "sex_classfemale", reporting_label] == "Female vs male")

qa <- rbindlist(checks)
fwrite(qa, file.path(output_root, "aggregate", "synthetic_checks.tsv"), sep = "\t")
writeLines(if (all(qa$pass)) "S405L_R1_SYNTHETIC_TESTS_PASS" else "S405L_R1_SYNTHETIC_TESTS_FAIL", file.path(output_root, "status.txt"))
if (!all(qa$pass)) stop("S405L_R1_SYNTHETIC_TEST_FAILURE", call. = FALSE)
