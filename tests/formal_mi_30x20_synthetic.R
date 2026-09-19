args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: OUTPUT_ROOT", call. = FALSE)
output_root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
dir.create(file.path(output_root, "aggregate"), recursive = TRUE, showWarnings = FALSE)

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][[1L]])
package_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
extra_library <- Sys.getenv("R9_EXTRA_R_LIB", unset = "")
if (nzchar(extra_library)) {
  if (!dir.exists(extra_library)) stop("EXTRA_R_LIBRARY_MISSING", call. = FALSE)
  .libPaths(c(extra_library, .libPaths()))
}
suppressPackageStartupMessages({
  library(data.table); library(mice); library(miceadds); library(posterior)
})
formal_root <- file.path(package_root, "code", "formal")
source(file.path(formal_root, "stage405_multilevel_core.R"))
source(file.path(formal_root, "stage405_formal_convergence_adapter.R"))
source(file.path(formal_root, "stage405_formal_imputation.R"))

checks <- list()
add_check <- function(check, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = check, pass = isTRUE(pass), detail = as.character(detail))
}

set.seed(20260907)
stays <- data.table(
  patient_key_internal = paste0("patient_", seq_len(30L)),
  stay_key_internal = paste0("stay_", seq_len(30L)),
  age_years_numeric = sample(30:88, 30L, replace = TRUE),
  age_interval = NA_character_,
  sex_class = sample(c("female", "male"), 30L, replace = TRUE),
  cci_value = sample(0:8, 30L, replace = TRUE)
)
stays[c(3L, 11L, 19L), age_years_numeric := NA_real_]
stays[c(6L, 17L), sex_class := NA_character_]
stays[c(8L, 22L), cci_value := NA_real_]
x <- stays[rep(seq_len(.N), each = 5L)]
x[, pair_uid := paste0("pair_", seq_len(.N))]
x[, `:=`(
  cohort = "MIMIC",
  paired_mean_saturation_c90 = rnorm(.N, 5, 2),
  sao2_anchor_icu_hours = runif(.N, 0, 72),
  sofa_total = pmax(0, pmin(24, round(rnorm(.N, 7, 3)))),
  pao2_candidate_value = pmax(30, rnorm(.N, 85, 15)),
  paco2_candidate_value = pmax(15, rnorm(.N, 40, 7)),
  ph_candidate_value = pmin(7.7, pmax(6.8, rnorm(.N, 7.36, 0.08))),
  lactate_candidate_value = pmax(0.2, rlnorm(.N, log(2), 0.4)),
  hb_candidate_value = pmax(3, rnorm(.N, 10.5, 1.5))
)]
for (field in c("sofa_total", "pao2_candidate_value", "paco2_candidate_value", "ph_candidate_value")) {
  x[sample.int(.N, 18L), (field) := NA_real_]
}
x[, bias_spo2_minus_sao2 := 0.2 * paired_mean_saturation_c90 + rnorm(.N)]
x <- s405m_derive_age_fields(x, "MIMIC")
eligible <- c(
  "age_value_per_10y", "age_topcoded_indicator", "sex_class", "cci_value",
  "sofa_total", "pao2_candidate_value", "paco2_candidate_value", "ph_candidate_value"
)

result <- s405f_run_imputation(x, "MIMIC", eligible, 30L, 20L, 20260907L, FALSE)
validation <- s405m_validate_imputation(x, result, "MIMIC")
domain <- s405m_imputed_domain_audit(result)

add_check("formal_result_exactly_30x20", result$m == 30L && result$maxit == 20L)
add_check("thirty_completed_datasets_created", length(result$completed) == 30L)
add_check("no_residual_missing", all(validation$variable$completed_no_missing))
add_check("observed_values_unchanged", all(validation$variable$observed_unchanged))
add_check("stay_values_constant", all(validation$stay_constancy$stay_constant))
add_check("patient_sex_constant", all(validation$patient_sex_constancy$patient_constant))
add_check("all_domains_pass", all(domain$domain_pass))
add_check("formal_adapter_applied", result$spec$formal_convergence_adapter_applied)
add_check("formal_gate_is_eligible", result$convergence$gate$formal_gate_eligible[[1L]])
add_check("trace_has_30_chains", all(result$convergence$completeness$chain_n == 30L))
add_check("trace_has_20_iterations", all(result$convergence$completeness$iteration_n == 20L))
add_check("no_lmm_object_created", !any(c("fit", "fits", "lmm", "pooled") %in% names(result)))
add_check("logged_events_exportable", is.data.table(result$logged_events))

wrong_design <- tryCatch(
  s405f_run_imputation(x, "MIMIC", eligible, 3L, 5L, 1L, FALSE),
  error = function(e) conditionMessage(e)
)
add_check("formal_wrapper_rejects_non30x20", identical(wrong_design, "S405F_FORMAL_DESIGN_MUST_BE_30X20"))

checks_dt <- rbindlist(checks)
fwrite(checks_dt, file.path(output_root, "aggregate", "synthetic_formal_checks.tsv"), sep = "\t")
fwrite(result$convergence$gate, file.path(output_root, "aggregate", "synthetic_convergence_gate.tsv"), sep = "\t")
fwrite(result$convergence$summary, file.path(output_root, "aggregate", "synthetic_convergence_summary.tsv"), sep = "\t")
writeLines(if (all(checks_dt$pass)) "S405F_SYNTHETIC_FORMAL_INTEGRATION_PASS" else "S405F_SYNTHETIC_FORMAL_INTEGRATION_FAIL", file.path(output_root, "status.txt"))
if (!all(checks_dt$pass)) stop("S405F_SYNTHETIC_FORMAL_INTEGRATION_FAILURE", call. = FALSE)
