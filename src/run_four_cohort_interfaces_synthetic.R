args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("USAGE_OUTPUT_DIRECTORY", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
candidate_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
output_root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
allowed_parent <- normalizePath(file.path(candidate_root, "runs"), winslash = "/", mustWork = FALSE)
if (!identical(dirname(output_root), allowed_parent)) stop("OUTPUT_OUTSIDE_CANDIDATE_RUNS", call. = FALSE)
if (dir.exists(output_root) || file.exists(output_root)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_root)) stop("OUTPUT_CREATE_FAILED", call. = FALSE)
writeLines("FOUR_COHORT_SYNTHETIC_RUNNING", file.path(output_root, "status.txt"))
options(error = function() {
  writeLines("FOUR_COHORT_SYNTHETIC_FAIL", file.path(output_root, "status.txt"))
  quit(save = "no", status = 1L)
})

module_root <- file.path(candidate_root, "code", "cores")
module <- data.frame(
  name = c("R04", "R05"),
  relative_path = c(
    "r04_event_core.R",
    "r05_pairing_core.R"
  ),
  expected_sha256 = c(
    "30D257FB34B83D137668B61C6AD5BDA3EB5B3682CC57AE0E350582EA7FC35ABD",
    "2143A8057D9A0DB7A62878CADF2DF103F2E43528301E4EA22DC0ADD210DFE812"
  ), stringsAsFactors = FALSE
)
paths <- file.path(module_root, module$relative_path)
if (any(!file.exists(paths))) stop("MODULE_SOURCE_MISSING", call. = FALSE)
module$observed_sha256 <- toupper(unname(tools::sha256sum(paths)))
if (any(module$observed_sha256 != module$expected_sha256)) stop("MODULE_SOURCE_HASH_MISMATCH", call. = FALSE)
for (path in paths) source(path, local = FALSE)

checks <- list()
add_check <- function(cohort, check, pass) {
  checks[[length(checks) + 1L]] <<- data.frame(cohort = cohort, check = check, pass = isTRUE(pass))
  if (!isTRUE(pass)) stop(paste0("CHECK_FAILED_", cohort, "_", check), call. = FALSE)
}
expect_error <- function(expr, code) {
  observed <- tryCatch({ force(expr); "NO_ERROR" }, error = function(e) conditionMessage(e))
  identical(observed, code)
}

# These are invented fixtures only. No patient or source-data files are opened.
cohorts <- c("Amsterdam", "eICU", "SICDB", "Lianyungang")
summary_rows <- list()
all_keys <- character()
for (cohort in cohorts) {
  raw_stays <- paste0("synthetic_", tolower(cohort), "_stay_", 1:2)
  raw_patients <- paste0("synthetic_", tolower(cohort), "_patient_", 1:2)
  bridge <- r04_build_p1_bridge(cohort, data.frame(stay_key = raw_stays, patient_key = raw_patients))
  add_check(cohort, "bridge_unique_and_cohort_scoped", !anyDuplicated(bridge$stay_key_internal) &&
    all(startsWith(bridge$stay_key_internal, paste0(tolower(cohort), "_stay_"))))
  all_keys <- c(all_keys, bridge$stay_key_internal)
  if (cohort == "Amsterdam") {
    normalized <- r04_amsterdam_normalize_sao2(c(12311, 8903, 11543), c(0.85, 90, 95))
    add_check(cohort, "fraction_normalized_only_for_item12311", identical(as.numeric(normalized$value), c(85, 90, 95)) &&
      identical(as.logical(normalized$fraction_normalized), c(TRUE, FALSE, FALSE)))
    add_check(cohort, "explicit_venous_comments_detected", identical(
      as.logical(r04_amsterdam_is_venous_comment(c("arterial", "venous sample", "veneus"))), c(FALSE, TRUE, TRUE)))
  }
  if (cohort == "eICU") {
    observationoffset <- c(100, 220, 340)
    labresultoffset <- observationoffset + c(1, 3, 30)
    add_check(cohort, "icu_relative_minute_interface", identical(as.numeric(labresultoffset - observationoffset), c(1, 3, 30)))
  }
  if (cohort == "SICDB") {
    icu_offset_sec <- 3600
    source_offset_sec <- icu_offset_sec + 60 * c(100, 220, 340)
    add_check(cohort, "seconds_to_icu_minutes", identical(as.numeric((source_offset_sec - icu_offset_sec) / 60), c(100, 220, 340)))
  }
  if (cohort == "Lianyungang") {
    add_check(cohort, "label_rules", identical(as.logical(r04_lyy_spo2_label_keep(c("SpO2", "血氧饱和度", "heart rate"))),
      c(TRUE, TRUE, FALSE)) && identical(as.logical(r04_lyy_sao2_label_keep(c("氧饱和度(SaO2)-动脉血", "氧饱和度(SaO2)-静脉血"))), c(TRUE, FALSE)))
    time_parsed <- r04_parse_lyy_utc_min(c("2026-01-01 00:00:00", "2026-01-01 00:01:00"))
    add_check(cohort, "source_time_and_value_parse", isTRUE(all.equal(diff(time_parsed), 1)) &&
      identical(as.numeric(r04_parse_lyy_value(c("85%", "90%"))), c(85, 90)))
  }
  make_staged <- function(modality) {
    x <- expand.grid(stay_index = 1:2, measurement = 1:3)
    stay_index <- x$stay_index
    j <- x$measurement
    base <- c(100, 220, 340)[j]
    is_spo <- modality == "SpO2"
    value <- c(85, 90, 95)[j] + if (is_spo) c(3, 1, -1)[j] else 0
    d <- data.frame(
      stay_raw = raw_stays[stay_index],
      stay_key_internal = bridge$stay_key_internal[stay_index],
      patient_key_internal = bridge$patient_key_internal[stay_index],
      event_time_min = as.numeric(base + if (is_spo) 0 else c(1, 3, 30)[j]),
      saturation_percent = as.numeric(value),
      source_family = paste0("synthetic_", tolower(cohort), "_", modality),
      source_priority = 1L, stringsAsFactors = FALSE
    )
    if (is_spo) {
      d <- rbind(d, d[1L, , drop = FALSE])
      extra <- d[1L, , drop = FALSE]
      extra$event_time_min <- 700
      extra$saturation_percent <- 69
      d <- rbind(d, extra)
    } else {
      extra <- d[1L, , drop = FALSE]
      extra$event_time_min <- 700
      extra$saturation_percent <- 101
      d <- rbind(d, extra)
    }
    d
  }
  artifact_assessable <- cohort != "SICDB"
  sp <- r04_finalize_events(make_staged("SpO2"), cohort, "SpO2",
    artifact_threshold = if (cohort == "Lianyungang") 18.286 else 11,
    artifact_assessable = artifact_assessable)
  sa <- r04_finalize_events(make_staged("SaO2"), cohort, "SaO2")
  r04_assert_controlled_event_output(sp$events)
  r04_assert_controlled_event_output(sa$events)
  add_check(cohort, "no_raw_fields", !any(c("stay_raw", "patient_raw") %in% names(sp$events)) &&
    !any(c("stay_raw", "patient_raw") %in% names(sa$events)))
  add_check(cohort, "duplicate_and_range", sp$aggregate$exact_duplicate_rows_removed[[1L]] == 1 &&
    sp$aggregate$below_70_events[[1L]] == 1 && sa$aggregate$above_100_events[[1L]] == 1)
  add_check(cohort, "artifact_mode", identical(sp$aggregate$artifact_qc_mode[[1L]],
    if (artifact_assessable) "P99_15_MIN_THREE_POINT" else "NOT_ASSESSABLE_TEMPORAL_RESOLUTION"))
  add_check(cohort, "invalid_staged_schema_rejected", expect_error(
    r04_finalize_events(data.frame(x = 1), cohort, "SpO2"), "R04_EVENT_STAGE_SCHEMA"))
  m60 <- r05_pair_cohort(as.data.frame(sp$events), as.data.frame(sa$events), cohort, 60)
  m5 <- r05_pair_cohort(as.data.frame(sp$events), as.data.frame(sa$events), cohort, 5)
  add_check(cohort, "independent_pair_windows", nrow(m60$final_pairs) == 6 && nrow(m5$final_pairs) == 4)
  add_check(cohort, "range_excluded", m60$aggregate$eligible_spo2_events[[1L]] == 6 &&
    m60$aggregate$eligible_sao2_events[[1L]] == 6)
  for (z in list(m60, m5)) {
    add_check(cohort, paste0("one_to_one_", z$aggregate$window_minutes[[1L]]),
      !anyDuplicated(paste(z$final_pairs$stay_key_internal, z$final_pairs$spo2_stable_event_ordinal)) &&
      !anyDuplicated(paste(z$final_pairs$stay_key_internal, z$final_pairs$sao2_stable_event_ordinal)))
  }
  shuffled <- r05_pair_cohort(as.data.frame(sp$events[nrow(sp$events):1, ]),
    as.data.frame(sa$events[nrow(sa$events):1, ]), cohort, 60)
  add_check(cohort, "pairing_deterministic", identical(r05_pair_signature(m60$final_pairs), r05_pair_signature(shuffled$final_pairs)))
  bad <- as.data.frame(sa$events)
  bad$patient_key_internal[[1L]] <- paste0(bad$patient_key_internal[[1L]], "_wrong")
  add_check(cohort, "cross_modality_patient_conflict_rejected", expect_error(
    r05_pair_cohort(as.data.frame(sp$events), bad, cohort, 60), "R05_STAY_PATIENT_INTERFACE_ERROR"))
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(cohort = cohort, synthetic_stays = 2L,
    eligible_spo2_events = m60$aggregate$eligible_spo2_events[[1L]],
    eligible_sao2_events = m60$aggregate$eligible_sao2_events[[1L]],
    pairs_60_min = nrow(m60$final_pairs), pairs_5_min = nrow(m5$final_pairs),
    artifact_qc_mode = sp$aggregate$artifact_qc_mode[[1L]], stringsAsFactors = FALSE)
}
add_check("ALL", "distinct_internal_stay_keys_between_cohorts", !anyDuplicated(all_keys))
check_table <- do.call(rbind, checks)
if (!all(check_table$pass)) stop("FOUR_COHORT_CHECKS_FAILED", call. = FALSE)
write.table(module, file.path(output_root, "source_receipt.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, summary_rows), file.path(output_root, "aggregate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(check_table, file.path(output_root, "checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines("FOUR_COHORT_SYNTHETIC_PASS_PENDING_INDEPENDENT_QA", file.path(output_root, "status.txt"))
