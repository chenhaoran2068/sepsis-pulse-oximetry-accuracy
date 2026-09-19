r211e_expected_cohorts <- function() c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")

r211e_protected_fields <- function() {
  c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid", "bias_spo2_minus_sao2", "spo2_percent", "sao2_percent")
}

r211e_require_columns <- function(x, required, code) {
  if (!all(required %in% names(x))) stop(code, call. = FALSE)
  invisible(TRUE)
}

r211e_validate_seed <- function(seed, code) {
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) || seed != as.integer(seed)) stop(code, call. = FALSE)
  invisible(TRUE)
}

r211e_validate_candidate <- function(candidate, cohort) {
  r211e_require_columns(candidate, r211e_protected_fields(), "R211E_CANDIDATE_INTERFACE_MISSING")
  if (nrow(candidate) == 0L || anyNA(candidate$cohort) || !all(as.character(candidate$cohort) == cohort)) stop("R211E_COHORT_IDENTITY_INVALID", call. = FALSE)
  if (anyNA(candidate$analysis_patient_key) || anyNA(candidate$analysis_stay_key) || anyNA(candidate$pair_uid) || anyDuplicated(candidate$pair_uid)) stop("R211E_PROTECTED_KEY_INTEGRITY_FAILURE", call. = FALSE)
  if (any(!is.finite(candidate$spo2_percent)) || any(!is.finite(candidate$sao2_percent)) || any(!is.finite(candidate$bias_spo2_minus_sao2))) stop("R211E_STRUCTURAL_OXYGEN_OR_BIAS_MISSING", call. = FALSE)
  if (!isTRUE(all.equal(candidate$bias_spo2_minus_sao2, candidate$spo2_percent - candidate$sao2_percent, check.attributes = FALSE))) stop("R211E_BIAS_DIRECTION_INVALID", call. = FALSE)
  invisible(TRUE)
}

r211e_assign_outer <- function(candidate, seed = 20260827L, development_fraction = 0.70) {
  r211e_validate_seed(seed, "R211E_OUTER_SEED_INVALID")
  if (length(development_fraction) != 1L || !is.finite(development_fraction) || development_fraction <= 0 || development_fraction >= 1) stop("R211E_OUTER_FRACTION_INVALID", call. = FALSE)
  patients <- sort(unique(as.character(candidate$analysis_patient_key)))
  if (length(patients) < 2L) stop("R211E_OUTER_PATIENT_COUNT_INSUFFICIENT", call. = FALSE)
  set.seed(as.integer(seed))
  ordered <- sample(patients, length(patients), replace = FALSE)
  n_development <- max(1L, min(length(patients) - 1L, as.integer(round(length(patients) * development_fraction))))
  data.frame(analysis_patient_key = ordered, outer_partition = c(rep("development", n_development), rep("test", length(patients) - n_development)), stringsAsFactors = FALSE)
}

r211e_assign_inner <- function(outer_allocation, seed = 20260828L, k = 5L) {
  r211e_require_columns(outer_allocation, c("analysis_patient_key", "outer_partition"), "R211E_OUTER_ALLOCATION_INTERFACE_MISSING")
  r211e_validate_seed(seed, "R211E_INNER_SEED_INVALID")
  if (length(k) != 1L || is.na(k) || k != 5L) stop("R211E_INNER_FOLD_COUNT_INVALID", call. = FALSE)
  development_patients <- sort(unique(as.character(outer_allocation$analysis_patient_key[outer_allocation$outer_partition == "development"])))
  if (length(development_patients) < k) stop("R211E_INNER_DEVELOPMENT_PATIENT_COUNT_INSUFFICIENT", call. = FALSE)
  set.seed(as.integer(seed))
  ordered <- sample(development_patients, length(development_patients), replace = FALSE)
  data.frame(analysis_patient_key = ordered, inner_fold = rep(seq_len(k), length.out = length(ordered)), stringsAsFactors = FALSE)
}

r211e_build_allocation <- function(candidate, cohort, outer_seed = 20260827L, inner_seed = 20260828L) {
  r211e_validate_candidate(candidate, cohort)
  outer <- r211e_assign_outer(candidate, outer_seed)
  inner <- r211e_assign_inner(outer, inner_seed)
  key_rows <- candidate[, c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid"), drop = FALSE]
  allocation <- merge(key_rows, outer, by = "analysis_patient_key", all.x = TRUE, sort = FALSE)
  allocation$inner_fold <- inner$inner_fold[match(as.character(allocation$analysis_patient_key), as.character(inner$analysis_patient_key))]
  allocation <- allocation[, c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid", "outer_partition", "inner_fold"), drop = FALSE]
  r211e_validate_allocation(candidate, allocation, cohort)
  allocation
}

r211e_validate_allocation <- function(candidate, allocation, cohort) {
  r211e_validate_candidate(candidate, cohort)
  r211e_require_columns(allocation, c("cohort", "analysis_patient_key", "analysis_stay_key", "pair_uid", "outer_partition", "inner_fold"), "R211E_ALLOCATION_INTERFACE_MISSING")
  if (nrow(allocation) != nrow(candidate) || anyDuplicated(allocation$pair_uid) || !setequal(as.character(allocation$pair_uid), as.character(candidate$pair_uid))) stop("R211E_PAIR_COVERAGE_INVALID", call. = FALSE)
  if (anyNA(allocation$cohort) || !all(as.character(allocation$cohort) == cohort)) stop("R211E_ALLOCATION_COHORT_INVALID", call. = FALSE)
  index <- match(as.character(candidate$pair_uid), as.character(allocation$pair_uid))
  if (anyNA(index) || !identical(as.character(candidate$analysis_patient_key), as.character(allocation$analysis_patient_key[index])) || !identical(as.character(candidate$analysis_stay_key), as.character(allocation$analysis_stay_key[index]))) stop("R211E_PROTECTED_KEY_LINKAGE_INVALID", call. = FALSE)
  patient_map <- unique(allocation[, c("analysis_patient_key", "outer_partition"), drop = FALSE])
  if (anyDuplicated(patient_map$analysis_patient_key) || any(!patient_map$outer_partition %in% c("development", "test")) || length(unique(patient_map$outer_partition)) != 2L) stop("R211E_OUTER_ALLOCATION_INVALID", call. = FALSE)
  if (any(is.na(allocation$inner_fold[allocation$outer_partition == "development"])) || any(!is.na(allocation$inner_fold[allocation$outer_partition == "test"]))) stop("R211E_INNER_TEST_SEPARATION_INVALID", call. = FALSE)
  inner_map <- unique(allocation[allocation$outer_partition == "development", c("analysis_patient_key", "inner_fold"), drop = FALSE])
  if (anyDuplicated(inner_map$analysis_patient_key) || any(!inner_map$inner_fold %in% seq_len(5L)) || !setequal(unique(inner_map$inner_fold), seq_len(5L))) stop("R211E_INNER_ALLOCATION_INVALID", call. = FALSE)
  invisible(TRUE)
}

r211e_allocation_summary <- function(candidate, allocation, cohort) {
  patient_partitions <- unique(allocation[, c("analysis_patient_key", "outer_partition"), drop = FALSE])
  stay_partitions <- unique(allocation[, c("analysis_stay_key", "outer_partition"), drop = FALSE])
  data.frame(
    cohort = cohort,
    input_pair_n = nrow(candidate),
    input_patient_n = length(unique(candidate$analysis_patient_key)),
    input_stay_n = length(unique(candidate$analysis_stay_key)),
    development_pair_n = sum(allocation$outer_partition == "development"),
    test_pair_n = sum(allocation$outer_partition == "test"),
    development_patient_n = sum(patient_partitions$outer_partition == "development"),
    test_patient_n = sum(patient_partitions$outer_partition == "test"),
    development_stay_n = sum(stay_partitions$outer_partition == "development"),
    test_stay_n = sum(stay_partitions$outer_partition == "test"),
    stringsAsFactors = FALSE
  )
}

r211e_inner_summary <- function(allocation, cohort) {
  development <- allocation[allocation$outer_partition == "development", , drop = FALSE]
  data.frame(
    cohort = cohort,
    inner_fold = seq_len(5L),
    patient_n = vapply(seq_len(5L), function(f) length(unique(development$analysis_patient_key[development$inner_fold == f])), integer(1)),
    pair_n = vapply(seq_len(5L), function(f) sum(development$inner_fold == f), integer(1)),
    stringsAsFactors = FALSE
  )
}

r211e_record_receipt <- function(role, path, parser = "base") {
  if (!file.exists(path)) stop("R211E_REQUIRED_INPUT_ABSENT", call. = FALSE)
  info <- file.info(path)
  data.frame(asset_role = role, path = normalizePath(path, winslash = "/"), bytes = as.numeric(info$size), checksum_algorithm = "sha256", checksum = digest::digest(file = path, algo = "sha256"), parser = parser, stringsAsFactors = FALSE)
}
