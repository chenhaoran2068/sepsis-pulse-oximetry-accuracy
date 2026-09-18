args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("USAGE_EVENT_RUN_AND_NEW_OUTPUT", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
run_parent <- normalizePath(file.path(root, "runs"), winslash = "/", mustWork = FALSE)
event_run <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (!identical(dirname(event_run), run_parent) || !identical(dirname(output), run_parent))
  stop("RUN_PATH_OUTSIDE_CANDIDATE", call. = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (!identical(readLines(file.path(event_run, "independent_qa_status.txt"), warn = FALSE),
               "STANDARDIZED_EVENT_INDEPENDENT_QA_PASS")) stop("EVENT_RUN_QA_REQUIRED", call. = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output)) stop("OUTPUT_CREATE_FAILED", call. = FALSE)
writeLines("CLINICAL_INPUT_DEMO_RUNNING", file.path(output, "status.txt"))
options(error = function() {
  writeLines("CLINICAL_INPUT_DEMO_FAIL", file.path(output, "status.txt"))
  quit(save = "no", status = 1L)
})

source(file.path(root, "code", "cores", "r04_event_core.R"), local = FALSE)
source(file.path(root, "code", "cores", "r05_pairing_core.R"), local = FALSE)
source(file.path(root, "code", "standardized_event_input_v1.R"), local = FALSE)
source(file.path(root, "code", "standardized_clinical_input_v1.R"), local = FALSE)
loaded <- si_load_input(file.path(event_run, "input"), "MIMIC")
pair60 <- r05_pair_cohort(loaded$spo2, loaded$sao2, "MIMIC", 60)$final_pairs
pair5 <- r05_pair_cohort(loaded$spo2, loaded$sao2, "MIMIC", 5)$final_pairs
stopifnot(nrow(pair60) == 80L, nrow(pair5) == 60L)

stays <- loaded$stays
stay_order <- match(stays$stay_key_internal, unique(stays$stay_key_internal))
stay_clinical <- data.frame(
  patient_key_internal = stays$patient_key_internal,
  stay_key_internal = stays$stay_key_internal,
  age_years_numeric = 30 + stay_order,
  age_interval = NA_character_,
  age_model_representation = "numeric",
  sex_class = ifelse(stay_order %% 2L, "female", "male"),
  cci_value = stay_order %% 6L,
  stringsAsFactors = FALSE
)
# The one invented patient with two stays has a consistent observed sex.
stay_clinical$sex_class[2L] <- stay_clinical$sex_class[1L]
metric_units <- c(pao2 = "mmHg", paco2 = "mmHg", ph = "pH", lactate = "mmol/L", hb = "g/dL")
metric_values <- c(pao2 = 78, paco2 = 43, ph = 7.34, lactate = 2.1, hb = 11.8)
make_clinical <- function(pairs, window) {
  x <- pairs[ci1_pair_keys()]
  x$pair_key_internal <- sprintf("invented_%s_pair_%04d", window, seq_len(nrow(pairs)))
  x$sofa_total <- as.numeric(seq_len(nrow(pairs)) %% 12L)
  x$sofa_status <- "available"
  for (name in names(metric_units)) {
    x[[paste0(name, "_candidate_value")]] <- metric_values[[name]]
    x[[paste0(name, "_canonical_unit")]] <- metric_units[[name]]
    x[[paste0(name, "_linkage_class")]] <- "same_specimen"
    x[[paste0(name, "_value_state")]] <- "eligible_pending_candidate_build"
  }
  for (name in ci1_treatment_names()) {
    x[[paste0(name, "_state")]] <- "unknown"
    x[[paste0(name, "_scope")]] <- "permitted_by_gate_r08"
  }
  x[ci1_pair_columns()]
}
clinical60 <- make_clinical(pair60, "60")
clinical5 <- make_clinical(pair5, "5")
input_dir <- file.path(output, "input")
dir.create(input_dir)
write_tsv <- function(x, name) utils::write.table(x, file.path(input_dir, name),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
write_tsv(stay_clinical, "stay_clinical.tsv")
write_tsv(clinical60, "pair_clinical_60.tsv")
write_tsv(clinical5, "pair_clinical_5.tsv")

checks <- list()
check <- function(name, value) {
  ok <- isTRUE(value)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = ok)
  if (!ok) stop(paste0("CHECK_FAILED_", name), call. = FALSE)
}
expect_error <- function(expr, code) identical(tryCatch({force(expr); "NO_ERROR"},
  error = function(e) conditionMessage(e)), code)
stay_read <- ci1_read_stays(file.path(input_dir, "stay_clinical.tsv"), loaded$stays)
clinical60_read <- ci1_read_pairs(file.path(input_dir, "pair_clinical_60.tsv"), pair60, stay_read, "MIMIC", 60)
clinical5_read <- ci1_read_pairs(file.path(input_dir, "pair_clinical_5.tsv"), pair5, stay_read, "MIMIC", 5)
analytic60 <- ci1_make_analytic_pairs(pair60, stay_read, clinical60_read, "MIMIC", 60)
analytic5 <- ci1_make_analytic_pairs(pair5, stay_read, clinical5_read, "MIMIC", 5)
check("stay_row_cardinality", nrow(stay_read) == 20L)
check("independent_pair_cardinality", nrow(analytic60) == 80L && nrow(analytic5) == 60L)
check("outcome_and_background", all(analytic60$bias_spo2_minus_sao2 ==
  analytic60$spo2_saturation_percent - analytic60$sao2_saturation_percent) &&
  all(analytic60$paired_mean_saturation_c90 ==
  (analytic60$spo2_saturation_percent + analytic60$sao2_saturation_percent) / 2 - 90))
check("stay_level_value_consistency", all(vapply(split(analytic60$cci_value, analytic60$stay_key_internal),
  function(x) length(unique(x)) == 1L, logical(1L))))
bad <- stay_read
bad$patient_key_internal[[1L]] <- "invented_wrong_patient"
check("wrong_stay_mapping_rejected", expect_error(ci1_read_stays(file.path(input_dir, "stay_clinical.tsv"), bad),
  "CI1_STAY_COVERAGE_MISMATCH"))
bad <- clinical60_read[-1L, , drop = FALSE]
check("missing_pair_rejected", expect_error(ci1_read_pairs(file.path(input_dir, "pair_clinical_60.tsv"),
  pair60[-1L, , drop = FALSE], stay_read, "MIMIC", 60), "CI1_PAIR_COVERAGE_MISMATCH"))
bad <- clinical60_read
bad$pair_key_internal[[2L]] <- bad$pair_key_internal[[1L]]
check("duplicate_pair_uid_rejected", expect_error(ci1_make_analytic_pairs(pair60, stay_read, bad, "MIMIC", 60),
  "CI1_ANALYTIC_CARDINALITY_INVALID"))
bad <- clinical60_read
bad$pao2_canonical_unit[[1L]] <- "kPa"
bad_file <- file.path(input_dir, "bad_unit.tsv")
write_tsv(bad, "bad_unit.tsv")
check("unit_mismatch_rejected", expect_error(ci1_read_pairs(bad_file, pair60, stay_read, "MIMIC", 60),
  "CI1_METRIC_UNIT_INVALID"))
bad <- clinical60_read
bad$ph_candidate_value[[1L]] <- NA_real_
write_tsv(bad, "bad_state.tsv")
check("state_value_mismatch_rejected", expect_error(ci1_read_pairs(file.path(input_dir, "bad_state.tsv"),
  pair60, stay_read, "MIMIC", 60), "CI1_METRIC_STATE_VALUE_MISMATCH"))
check("missing_input_rejected", expect_error(ci1_read_table(file.path(input_dir, "absent.tsv"),
  ci1_pair_columns()), "CI1_INPUT_FILE_ABSENT"))
check("window_mismatch_rejected", expect_error(ci1_read_pairs(file.path(input_dir, "pair_clinical_60.tsv"),
  pair60, stay_read, "MIMIC", 5), "CI1_MATCHED_PAIR_SCOPE_INVALID"))

summary <- data.frame(window_minutes = c(60L, 5L), stay_n = nrow(stay_read),
  patient_n = length(unique(stay_read$patient_key_internal)),
  analytic_pair_n = c(nrow(analytic60), nrow(analytic5)))
utils::write.table(summary, file.path(output, "aggregate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(do.call(rbind, checks), file.path(output, "checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("CLINICAL_INPUT_DEMO_PASS", file.path(output, "status.txt"))
