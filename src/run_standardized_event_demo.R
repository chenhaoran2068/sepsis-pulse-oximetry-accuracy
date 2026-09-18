args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("USAGE_NEW_OUTPUT_DIRECTORY", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
candidate_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
output_root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
allowed_parent <- normalizePath(file.path(candidate_root, "runs"), winslash = "/", mustWork = FALSE)
if (!identical(dirname(output_root), allowed_parent)) stop("OUTPUT_OUTSIDE_CANDIDATE_RUNS", call. = FALSE)
if (file.exists(output_root) || dir.exists(output_root)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_root)) stop("OUTPUT_CREATE_FAILED", call. = FALSE)
writeLines("STANDARDIZED_EVENT_DEMO_RUNNING", file.path(output_root, "status.txt"))
options(error = function() {
  writeLines("STANDARDIZED_EVENT_DEMO_FAIL", file.path(output_root, "status.txt"))
  quit(save = "no", status = 1L)
})

core_root <- file.path(candidate_root, "code", "cores")
core_paths <- file.path(core_root, c("r04_event_core.R", "r05_pairing_core.R"))
expected <- c("30D257FB34B83D137668B61C6AD5BDA3EB5B3682CC57AE0E350582EA7FC35ABD",
  "2143A8057D9A0DB7A62878CADF2DF103F2E43528301E4EA22DC0ADD210DFE812")
if (any(!file.exists(core_paths)) || any(toupper(unname(tools::sha256sum(core_paths))) != expected)) stop("FORMAL_CORE_HASH_MISMATCH", call. = FALSE)
for (path in core_paths) source(path, local = FALSE)
source(file.path(candidate_root, "code", "standardized_event_input_v1.R"), local = FALSE)
suppressPackageStartupMessages(library(data.table))

check_rows <- list()
check <- function(name, condition) {
  pass <- isTRUE(condition)
  check_rows[[length(check_rows) + 1L]] <<- data.frame(check = name, pass = pass)
  if (!pass) stop(paste0("CHECK_FAILED_", name), call. = FALSE)
}
expect_error <- function(expr, code) {
  observed <- tryCatch({ force(expr); "NO_ERROR" }, error = function(e) conditionMessage(e))
  identical(observed, code)
}

# Invented source keys and measurements; none is sampled from a study cohort.
n_stay <- 20L
raw_stay <- sprintf("invented_stay_%03d", seq_len(n_stay))
raw_patient <- sprintf("invented_patient_%03d", c(1L, 1L, seq.int(3L, n_stay)))
bridge <- r04_build_p1_bridge("MIMIC", data.frame(stay_key = raw_stay, patient_key = raw_patient))
stays <- bridge[, c("stay_key_internal", "patient_key_internal")]
make_staged <- function(modality) {
  blocks <- vector("list", n_stay * 4L)
  k <- 0L
  for (i in seq_len(n_stay)) for (j in seq_len(4L)) {
    k <- k + 1L
    sao <- c(86, 89, 92, 96)[[j]]
    spo <- sao + c(3, 1, -1, 0)[[j]]
    blocks[[k]] <- data.frame(stay_raw = raw_stay[[i]],
      stay_key_internal = bridge$stay_key_internal[[i]],
      patient_key_internal = bridge$patient_key_internal[[i]],
      event_time_min = 100 + (j - 1L) * 120 + if (modality == "SpO2") c(1, 4, 30, -2)[[j]] else 0,
      saturation_percent = if (modality == "SpO2") spo else sao,
      source_family = paste0("invented_", modality), source_priority = 1L)
  }
  out <- do.call(rbind, blocks)
  extra <- out[1L, , drop = FALSE]
  extra$event_time_min <- 700
  extra$saturation_percent <- if (modality == "SpO2") 69 else 101
  rbind(out, extra)
}
spo2 <- as.data.frame(r04_finalize_events(make_staged("SpO2"), "MIMIC", "SpO2", artifact_assessable = FALSE)$events)
sao2 <- as.data.frame(r04_finalize_events(make_staged("SaO2"), "MIMIC", "SaO2")$events)
input_dir <- file.path(output_root, "input")
dir.create(input_dir, recursive = FALSE)
write_tsv <- function(data, name) utils::write.table(data, file.path(input_dir, name), sep = "\t",
  row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
write_tsv(stays, "stays.tsv")
write_tsv(spo2, "spo2_events.tsv")
write_tsv(sao2, "sao2_events.tsv")

loaded <- si_load_input(input_dir, "MIMIC")
check("synthetic_stay_count", nrow(loaded$stays) == n_stay)
check("same_patient_multiple_stays_allowed", length(unique(loaded$stays$patient_key_internal)) == n_stay - 1L)
check("source_event_count", nrow(loaded$spo2) == 4L * n_stay + 1L && nrow(loaded$sao2) == 4L * n_stay + 1L)
check("range_flags_preserved", sum(!loaded$spo2$analysis_eligible_range) == 1L && sum(!loaded$sao2$analysis_eligible_range) == 1L)
paired60 <- r05_pair_cohort(loaded$spo2, loaded$sao2, "MIMIC", 60)
paired5 <- r05_pair_cohort(loaded$spo2, loaded$sao2, "MIMIC", 5)
check("independent_pairing_windows", nrow(paired60$final_pairs) == 4L * n_stay && nrow(paired5$final_pairs) == 3L * n_stay)
check("patient_stay_mapping", all(paired60$final_pairs$patient_key_internal %in% stays$patient_key_internal))

raw_spo <- si_read_tsv(file.path(input_dir, "spo2_events.tsv"), si_required_event_columns, "SPO2")
bad <- raw_spo
bad$stay_key_internal[[1L]] <- "unknown_stay"
check("unknown_stay_rejected", expect_error(si_validate_events(bad, loaded$stays, "MIMIC", "SpO2"), "SI_EVENT_STAY_NOT_IN_ROSTER"))
bad <- raw_spo
bad$patient_key_internal[[1L]] <- "wrong_patient"
check("patient_mismatch_rejected", expect_error(si_validate_events(bad, loaded$stays, "MIMIC", "SpO2"), "SI_EVENT_PATIENT_MISMATCH"))
bad <- raw_spo
bad$range_status[[1L]] <- "below_70"
check("range_mismatch_rejected", expect_error(si_validate_events(bad, loaded$stays, "MIMIC", "SpO2"), "SI_RANGE_VALUE_MISMATCH"))
bad <- raw_spo
bad$event_time_min[[1L]] <- "later"
check("invalid_time_rejected", expect_error(si_validate_events(bad, loaded$stays, "MIMIC", "SpO2"), "SI_EVENT_TIME_INVALID"))
bad <- raw_spo
bad$analysis_eligible_range[[1L]] <- "yes"
check("invalid_boolean_rejected", expect_error(si_validate_events(bad, loaded$stays, "MIMIC", "SpO2"), "SI_ELIGIBLE_FLAG_INVALID"))
bad <- raw_spo
bad$stable_event_ordinal[[2L]] <- bad$stable_event_ordinal[[1L]]
check("duplicate_event_ordinal_rejected", expect_error(si_validate_events(bad, loaded$stays, "MIMIC", "SpO2"), "R05_EVENT_ORDINAL_DUPLICATE"))
bad_stays <- rbind(loaded$stays, loaded$stays[1L, ])
check("duplicate_stay_rejected", expect_error(si_validate_stays(bad_stays), "SI_STAY_KEY_DUPLICATE"))
check("missing_input_rejected", expect_error(si_read_tsv(file.path(input_dir, "absent.tsv"), si_required_stay_columns, "STAYS"), "SI_STAYS_FILE_MISSING"))
check("empty_roster_rejected", expect_error(si_validate_stays(loaded$stays[0, ]), "SI_STAYS_EMPTY"))

aggregate <- data.frame(window_minutes = c(60L, 5L), stay_n = n_stay,
  patient_n = length(unique(stays$patient_key_internal)),
  final_pair_n = c(nrow(paired60$final_pairs), nrow(paired5$final_pairs)),
  eligible_spo2_n = c(paired60$aggregate$eligible_spo2_events, paired5$aggregate$eligible_spo2_events),
  eligible_sao2_n = c(paired60$aggregate$eligible_sao2_events, paired5$aggregate$eligible_sao2_events))
utils::write.table(aggregate, file.path(output_root, "aggregate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(do.call(rbind, check_rows), file.path(output_root, "checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("STANDARDIZED_EVENT_DEMO_PASS", file.path(output_root, "status.txt"))
cat("STANDARDIZED_EVENT_DEMO_PASS\n")
