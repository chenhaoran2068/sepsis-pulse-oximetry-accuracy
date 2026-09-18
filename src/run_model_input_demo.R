args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("USAGE_EVENT_RUN_CLINICAL_RUN_NEW_OUTPUT", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
run_parent <- normalizePath(file.path(root, "runs"), winslash = "/", mustWork = FALSE)
event_run <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
clinical_run <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
if (!identical(dirname(event_run), run_parent) ||
    !identical(dirname(clinical_run), run_parent) ||
    !identical(dirname(output), run_parent)) stop("RUN_PATH_OUTSIDE_CANDIDATE", call. = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (!identical(readLines(file.path(event_run, "independent_qa_status.txt"), warn = FALSE),
               "STANDARDIZED_EVENT_INDEPENDENT_QA_PASS") ||
    !identical(readLines(file.path(clinical_run, "independent_qa_status.txt"), warn = FALSE),
               "CLINICAL_INPUT_INDEPENDENT_QA_PASS")) stop("UPSTREAM_QA_REQUIRED", call. = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output)) stop("OUTPUT_CREATE_FAILED", call. = FALSE)
writeLines("MODEL_INPUT_DEMO_RUNNING", file.path(output, "status.txt"))
options(error = function() {
  writeLines("MODEL_INPUT_DEMO_FAIL", file.path(output, "status.txt"))
  quit(save = "no", status = 1L)
})

for (code in c("r04_event_core.R", "r05_pairing_core.R", "r211b_core.R",
               "stage405_mi_demo_core.R"))
  source(file.path(root, "code", "cores", code), local = FALSE)
source(file.path(root, "code", "standardized_event_input_v1.R"), local = FALSE)
source(file.path(root, "code", "standardized_clinical_input_v1.R"), local = FALSE)
source(file.path(root, "code", "model_input_adapter_v1.R"), local = FALSE)
loaded <- si_load_input(file.path(event_run, "input"), "MIMIC")
stays <- ci1_read_stays(file.path(clinical_run, "input", "stay_clinical.tsv"), loaded$stays)
checks <- list()
check <- function(name, value) {
  ok <- isTRUE(value)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = ok)
  if (!ok) stop(paste0("CHECK_FAILED_", name), call. = FALSE)
}
expect_error <- function(expr, code) identical(tryCatch({force(expr); "NO_ERROR"},
  error = function(e) conditionMessage(e)), code)
rows <- list()
for (window in c(60L, 5L)) {
  pair <- r05_pair_cohort(loaded$spo2, loaded$sao2, "MIMIC", window)$final_pairs
  clinical <- ci1_read_pairs(file.path(clinical_run, "input",
    paste0("pair_clinical_", window, ".tsv")), pair, stays, "MIMIC", window)
  analytic <- ci1_make_analytic_pairs(pair, stays, clinical, "MIMIC", window)
  lmm <- mi1_lmm_frame(analytic, pair, "MIMIC", window)
  check(paste0("lmm_anchor_", window), isTRUE(all.equal(lmm$sao2_anchor_icu_hours,
    pair$sao2_relative_icu_minutes / 60, check.attributes = FALSE)))
  check(paste0("lmm_structural_", window), isTRUE(all.equal(
    lmm$paired_mean_saturation_c90,
    (analytic$spo2_saturation_percent + analytic$sao2_saturation_percent) / 2 - 90,
    check.attributes = FALSE)))
  check(paste0("lmm_stay_age_", window),
    all(lmm$age_value_per_10y == analytic$age_years_numeric / 10))
  ledger <- s405m_coverage_ledger(lmm, "MIMIC")
  check(paste0("coverage_ledger_", window),
    nrow(ledger) == length(c(s405m_stay_fields("MIMIC"), s405m_pair_fields())))
  rows[[length(rows) + 1L]] <- data.frame(window_minutes = window,
    pair_n = nrow(analytic), patient_n = length(unique(analytic$patient_key_internal)),
    lmm_field_n = ncol(lmm))
  if (window == 60L) {
    rf <- mi1_rf_candidate(analytic, "MIMIC")
    check("rf_pair_identity", identical(as.character(rf$pair_uid),
      as.character(analytic$pair_key_internal)))
    check("rf_no_time_or_source_predictors",
      !any(grepl("relative_icu|lag|direction|artifact|source_label|source_code", names(rf))))
    check("rf_formula_identity", all(rf$bias_spo2_minus_sao2 ==
      rf$spo2_percent - rf$sao2_percent))
    check("rf_5_minute_rejected", expect_error(mi1_rf_candidate(
      within(analytic, pair_window_minutes <- 5L), "MIMIC"),
      "MI1_RF_REQUIRES_60_MINUTE_PAIRS"))
  } else {
    check("rf_sensitivity_window_rejected", expect_error(mi1_rf_candidate(
      analytic, "MIMIC"), "MI1_RF_REQUIRES_60_MINUTE_PAIRS"))
  }
}
summary <- do.call(rbind, rows)
utils::write.table(summary, file.path(output, "aggregate.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE)
utils::write.table(do.call(rbind, checks), file.path(output, "checks.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE)
writeLines("MODEL_INPUT_DEMO_PASS", file.path(output, "status.txt"))
