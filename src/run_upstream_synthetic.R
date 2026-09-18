args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) stop("USAGE_OUTPUT_DIRECTORY_AND_OPTIONAL_PATIENT_COUNT", call. = FALSE)
n_patient <- if (length(args) == 2L) suppressWarnings(as.integer(args[[2L]])) else 60L
if (is.na(n_patient) || n_patient < 20L || n_patient > 500L) stop("SYNTHETIC_PATIENT_COUNT_INVALID", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("SCRIPT_PATH_UNAVAILABLE", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
candidate_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
output_root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
if (dir.exists(output_root) || file.exists(output_root)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
allowed_parent <- normalizePath(file.path(candidate_root, "runs"), winslash = "/", mustWork = FALSE)
if (!identical(dirname(output_root), allowed_parent)) stop("OUTPUT_OUTSIDE_CANDIDATE_RUNS", call. = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_root)) stop("OUTPUT_CREATE_FAILED", call. = FALSE)
writeLines("SYNTHETIC_UPSTREAM_RUNNING", file.path(output_root, "status.txt"))
options(error = function() {
  writeLines("SYNTHETIC_UPSTREAM_FAIL", file.path(output_root, "status.txt"))
  quit(save = "no", status = 1L)
})

module_spec <- data.frame(
  module = c("R04", "R05", "Stage403", "Stage404", "Stage404B", "Stage404C"),
  relative_path = c(
    "r04_event_core.R",
    "r05_pairing_core.R",
    "agreement_core.R",
    "stage404_core.R",
    "stage404b_core.R",
    "stage404c_core.R"
  ),
  expected_sha256 = c(
    "30D257FB34B83D137668B61C6AD5BDA3EB5B3682CC57AE0E350582EA7FC35ABD",
    "2143A8057D9A0DB7A62878CADF2DF103F2E43528301E4EA22DC0ADD210DFE812",
    "1CAB7C999D74EF93AECD45E83EA6FE3D1FFCDC2FCA482345A5682CF2825782BF",
    "3BE2701C3DFDC213CDC2E3ADB445FA6148A179CCD5CE4D2785F694CCBC3DC726",
    "1434274EC5B6B018672BABBB1BDDA1A27052531BCED25BD63BB56C24599DC82D",
    "EA577D4446557F0C57F08A47B98F60AC3CF7AC3535E73DF4ED48CBF6AD3E673B"
  ), stringsAsFactors = FALSE
)
module_root <- file.path(candidate_root, "code", "cores")
module_paths <- file.path(module_root, module_spec$relative_path)
if (any(!file.exists(module_paths))) stop("MODULE_SOURCE_MISSING", call. = FALSE)
module_spec$observed_sha256 <- toupper(unname(tools::sha256sum(module_paths)))
if (any(module_spec$observed_sha256 != module_spec$expected_sha256)) stop("MODULE_SOURCE_HASH_MISMATCH", call. = FALSE)
for (path in module_paths) source(path, local = FALSE)
suppressPackageStartupMessages(library(data.table))

checks <- list()
add_check <- function(name, pass) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(pass), stringsAsFactors = FALSE)
  if (!isTRUE(pass)) stop(paste0("CHECK_FAILED_", name), call. = FALSE)
}
expect_error <- function(expr, code) {
  observed <- tryCatch({ force(expr); "NO_ERROR" }, error = function(e) conditionMessage(e))
  identical(observed, code)
}

# All keys and measurements below are invented in this script, not sampled from a cohort.
stay_raw <- sprintf("synthetic_stay_%03d", seq_len(n_patient))
patient_raw <- sprintf("synthetic_patient_%03d", seq_len(n_patient))
bridge <- r04_build_p1_bridge("MIMIC", data.frame(stay_key = stay_raw, patient_key = patient_raw))
make_staged <- function(modality) {
  rows <- vector("list", n_patient * 4L)
  k <- 0L
  for (i in seq_len(n_patient)) for (j in seq_len(4L)) {
    k <- k + 1L
    sao <- c(86, 89, 92, 96)[[j]] + ((i %% 3L) - 1L) * 0.2
    delta <- c(3, 1, -1, 0)[[j]] + ((i %% 7L) - 3L) * 0.6 +
      1.2 * sin(i * 0.47) + 0.15 * sin(i * j * 0.73) +
      if (j == 1L && i %% 5L == 0L) 4 else 0
    sat <- if (modality == "SpO2") sao + delta else sao
    minute <- 100 + (j - 1L) * 120 + ((i * 13L + j * j * 3L) %% 19L) +
      if (modality == "SpO2") c(1, 4, 30, -2)[[j]] else 0
    rows[[k]] <- data.frame(stay_raw = stay_raw[[i]], stay_key_internal = bridge$stay_key_internal[[i]],
      patient_key_internal = bridge$patient_key_internal[[i]], event_time_min = as.numeric(minute),
      saturation_percent = as.numeric(sat), source_family = paste0("synthetic_", modality),
      source_priority = 1L, stringsAsFactors = FALSE)
  }
  x <- do.call(rbind, rows)
  if (modality == "SpO2") {
    x <- rbind(x, x[1L, , drop = FALSE]) # exact duplicate
    extra <- x[1L, , drop = FALSE]
    extra$event_time_min <- 700
    extra$saturation_percent <- 69
    x <- rbind(x, extra)
  } else {
    extra <- x[1L, , drop = FALSE]
    extra$event_time_min <- 700
    extra$saturation_percent <- 101
    x <- rbind(x, extra)
  }
  x
}
sp <- r04_finalize_events(make_staged("SpO2"), "MIMIC", "SpO2", artifact_assessable = FALSE)
sa <- r04_finalize_events(make_staged("SaO2"), "MIMIC", "SaO2")
r04_assert_controlled_event_output(sp$events)
r04_assert_controlled_event_output(sa$events)
add_check("r04_duplicate_removed", sp$aggregate$exact_duplicate_rows_removed[[1L]] == 1)
add_check("r04_range_status", sp$aggregate$below_70_events[[1L]] == 1 && sa$aggregate$above_100_events[[1L]] == 1)
add_check("r04_no_raw_fields", !any(c("stay_raw", "patient_raw") %in% names(sp$events)))
add_check("r04_invalid_schema_rejected", expect_error(r04_finalize_events(data.frame(x = 1), "MIMIC", "SpO2"), "R04_EVENT_STAGE_SCHEMA"))

m60 <- r05_pair_cohort(as.data.frame(sp$events), as.data.frame(sa$events), "MIMIC", 60)
m5 <- r05_pair_cohort(as.data.frame(sp$events), as.data.frame(sa$events), "MIMIC", 5)
add_check("r05_independent_window_counts", nrow(m60$final_pairs) == 4L * n_patient && nrow(m5$final_pairs) == 3L * n_patient)
add_check("r05_range_excluded", m60$aggregate$eligible_spo2_events[[1L]] == 4L * n_patient && m60$aggregate$eligible_sao2_events[[1L]] == 4L * n_patient)
for (z in list(m60, m5)) {
  add_check(paste0("r05_one_to_one_", z$aggregate$window_minutes[[1L]]),
    !anyDuplicated(paste(z$final_pairs$stay_key_internal, z$final_pairs$spo2_stable_event_ordinal)) &&
    !anyDuplicated(paste(z$final_pairs$stay_key_internal, z$final_pairs$sao2_stable_event_ordinal)))
}
add_check("r05_maximum_cardinality", m60$aggregate$final_retained_pair_count[[1L]] == 4L * n_patient)
add_check("r05_minimum_total_lag", abs(sum(m60$final_pairs$absolute_lag_minutes) - n_patient * 37) < 1e-8)
shuffled <- r05_pair_cohort(as.data.frame(sp$events[sample.int(nrow(sp$events)), ]),
  as.data.frame(sa$events[sample.int(nrow(sa$events)), ]), "MIMIC", 60)
add_check("r05_shuffle_deterministic", identical(r05_pair_signature(m60$final_pairs), r05_pair_signature(shuffled$final_pairs)))
add_check("r05_missing_column_rejected", expect_error(r05_pair_cohort(as.data.frame(sp$events)[, setdiff(names(sp$events), "event_time_min")],
  as.data.frame(sa$events), "MIMIC", 60), "R05_EVENT_SCHEMA_ERROR"))

# A separate tiny case tests the tertiary stable-ordinal tie without touching real records.
small_event <- function(time, value, ordinal, modality) data.frame(
  stay_key_internal = "synthetic_tie_stay", patient_key_internal = "synthetic_tie_patient",
  event_time_min = as.numeric(time), saturation_percent = as.numeric(value),
  source_family = paste0("synthetic_", modality), source_priority = 1L,
  range_status = "eligible_70_100", analysis_eligible_range = TRUE,
  artifact_qc_status = "assessed_not_flagged", possible_transient_artifact = FALSE,
  stable_event_ordinal = as.numeric(ordinal), stringsAsFactors = FALSE)
tie_sp <- rbind(small_event(-1, 91, 2, "SpO2"), small_event(1, 92, 1, "SpO2"))
tie_sa <- small_event(0, 90, 1, "SaO2")
tie <- r05_pair_cohort(tie_sp, tie_sa, "MIMIC", 5)
add_check("r05_tertiary_tiebreak", nrow(tie$final_pairs) == 1L && tie$final_pairs$spo2_stable_event_ordinal[[1L]] == 1L)
artifact_sp <- rbind(small_event(0, 90, 1, "SpO2"), small_event(3, 100, 2, "SpO2"), small_event(6, 90, 3, "SpO2"))
artifact_sp$stay_key_internal <- "synthetic_artifact_stay"
artifact_sp$patient_key_internal <- "synthetic_artifact_patient"
artifact_sp$possible_transient_artifact[[2L]] <- TRUE
artifact_sp$artifact_qc_status[[2L]] <- "possible_transient_artifact"
artifact_sa <- small_event(3, 96, 1, "SaO2")
artifact_sa$stay_key_internal <- "synthetic_artifact_stay"
artifact_sa$patient_key_internal <- "synthetic_artifact_patient"
artifact <- r05_pair_cohort(artifact_sp, artifact_sa, "MIMIC", 5)
add_check("r05_artifact_removed_no_rematch", nrow(artifact$pre_artifact_pairs) == 1L && nrow(artifact$final_pairs) == 0L)

to_downstream <- function(pairs, window) {
  x <- data.table::as.data.table(pairs)
  x[, pair_key_internal := sprintf("synthetic_pair_%s_%04d", window, seq_len(.N))]
  x
}
run_downstream <- function(pairs, window, seed) {
  x <- to_downstream(pairs, window)
  validate_pair_input(x, "MIMIC", window)
  s404_validate_full_pairs(x, "MIMIC", window)
  z <- prepare_agreement_data(x)
  fit <- summarize_overall(z, b = 40L, seed = seed)
  add_check(paste0("stage403_finite_", window), all(is.finite(c(fit$fit$mean_bias, fit$fit$lower_loa, fit$fit$upper_loa, fit$dist$pair_arms))))
  add_check(paste0("stage403_bias_identity_", window), abs(fit$dist$pair_mean - mean(x$spo2_saturation_percent - x$sao2_saturation_percent)) < 1e-10)
  add_check(paste0("stage403_bootstrap_", window), sum(fit$boot$success) >= 38L)
  forward88 <- s404_overall_pattern(x, 88, 100L, seed + 1L)
  forward92 <- s404_overall_pattern(x, 92, 100L, seed + 2L)
  reverse88 <- s404b_reverse_pattern(x, 88, 100L, seed + 3L)
  reverse92 <- s404b_reverse_pattern(x, 92, 100L, seed + 4L)
  patient88 <- s404c_pattern_summary(x, 88)
  add_check(paste0("stage404_subset_", window), forward92$pair_numerator_n[[1L]] <= forward88$pair_numerator_n[[1L]])
  add_check(paste0("stage404_reverse_subset_", window), reverse92$conditional_numerator_sao2_lt88_n[[1L]] <= reverse88$conditional_numerator_sao2_lt88_n[[1L]])
  add_check(paste0("stage404_patient_once_", window), patient88$patient_summary$affected_patient_n[[1L]] <= n_patient)
  add_check(paste0("stage404_threshold_denominators_", window), forward88$pair_denominator_sao2_lt88_n[[1L]] == n_patient && reverse88$conditional_denominator_spo2_ge_n[[1L]] <= nrow(x))
  day <- data.table::copy(x)[sao2_relative_icu_minutes >= 0 & sao2_relative_icu_minutes <= 10080]
  day[, `:=`(relative_icu_day = pmin(7L, as.integer(floor(sao2_relative_icu_minutes / 1440) + 1L)),
    spo2_display_stratum = s404_source_stratum_code(spo2_saturation_percent),
    sao2_lt88 = sao2_saturation_percent < 88)]
  s404_validate_q21_pairs(day, "MIMIC", window)
  cell <- s404_day_stratum_cell(day[spo2_display_stratum == "70%-<88%"], 40L, seed + 5L)
  add_check(paste0("stage404_day_cell_", window), cell$numerator_n[[1L]] <= cell$denominator_n[[1L]])
  add_check(paste0("stage404_sparse_cell_", window), s404_day_stratum_cell(day[0], 40L, seed + 6L)$cell_status[[1L]] == "NO_DENOMINATOR")
  data.frame(window = window, final_pair_n = nrow(x), patient_n = uniqueN(x$patient_key_internal),
    mean_bias = fit$fit$mean_bias, lower_loa = fit$fit$lower_loa, upper_loa = fit$fit$upper_loa,
    arms = fit$dist$pair_arms, low_sao2_n = forward88$pair_denominator_sao2_lt88_n[[1L]],
    occult_n = forward88$pair_numerator_n[[1L]], severe_occult_n = forward92$pair_numerator_n[[1L]],
    reverse88_denominator_n = reverse88$conditional_denominator_spo2_ge_n[[1L]],
    reverse88_numerator_n = reverse88$conditional_numerator_sao2_lt88_n[[1L]],
    affected_patient_n = patient88$patient_summary$affected_patient_n[[1L]], stringsAsFactors = FALSE)
}
set.seed(20260916L)
aggregate <- rbind(run_downstream(m60$final_pairs, "M60", 11L), run_downstream(m5$final_pairs, "M5", 31L))
add_check("all_checks_pass", all(vapply(checks, function(z) z$pass[[1L]], logical(1))))

write.table(module_spec, file.path(output_root, "source_receipt.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(synthetic_patient_n = n_patient, synthetic_pairs_per_patient_m60 = 4L,
  synthetic_pairs_per_patient_m5 = 3L), file.path(output_root, "fixture_spec.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, checks), file.path(output_root, "checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(aggregate, file.path(output_root, "aggregate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(list(M60 = m60$final_pairs, M5 = m5$final_pairs), file.path(output_root, "synthetic_pairs_INTERNAL_ONLY.rds"))
writeLines("SYNTHETIC_UPSTREAM_PASS", file.path(output_root, "status.txt"))
cat("SYNTHETIC_UPSTREAM_PASS\n")
