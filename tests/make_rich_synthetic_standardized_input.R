args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: INVENTED_PAIR_RDS NEW_OUTPUT_DIR", call. = FALSE)
pair_source <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)
source(file.path(root, "code", "cores", "r05_pairing_core.R"))
source(file.path(root, "code", "standardized_event_input_v1.R"))
source(file.path(root, "code", "standardized_clinical_input_v1.R"))
source(file.path(root, "code", "cores", "r04_event_core.R"))
source(file.path(root, "code", "cores", "r211b_core.R"))
source(file.path(root, "code", "formal", "stage405_multilevel_core.R"))
source(file.path(root, "code", "model_input_adapter_v1.R"))
pairs <- readRDS(pair_source)
if (!is.list(pairs) || !all(c("M60", "M5") %in% names(pairs)) ||
    nrow(pairs$M60) < 800L || nrow(pairs$M5) != 0.75 * nrow(pairs$M60) ||
    !all(grepl("^mimic_", as.character(pairs$M60$patient_key_internal))))
  stop("INVENTED_PAIR_SOURCE_INVALID", call. = FALSE)
full <- pairs$M60
stays <- unique(as.data.frame(full[, c("stay_key_internal", "patient_key_internal")]))
make_events <- function(modality) {
  p <- paste0(modality, "_")
  result <- data.frame(
    stay_key_internal = full$stay_key_internal,
    patient_key_internal = full$patient_key_internal,
    event_time_min = full[[paste0(p, "relative_icu_minutes")]],
    saturation_percent = full[[paste0(p, "saturation_percent")]],
    source_family = paste0("invented_", modality),
    source_priority = 1L,
    range_status = "eligible_70_100",
    analysis_eligible_range = TRUE,
    artifact_qc_status = if (modality == "spo2") "not_assessable_temporal_resolution" else
      "not_applicable_non_spo2",
    possible_transient_artifact = FALSE,
    stable_event_ordinal = full[[paste0(p, "stable_event_ordinal")]],
    stringsAsFactors = FALSE)
  result
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) utils::write.table(x, file.path(output, name), sep = "\t",
                                                   row.names = FALSE, quote = FALSE, na = "")
write_tsv(stays, "stays.tsv")
write_tsv(make_events("spo2"), "spo2_events.tsv")
write_tsv(make_events("sao2"), "sao2_events.tsv")
loaded <- si_load_input(output, "MIMIC")
stay_i <- seq_len(nrow(stays))
stay_clinical <- data.frame(
  patient_key_internal = stays$patient_key_internal,
  stay_key_internal = stays$stay_key_internal,
  age_years_numeric = 31 + (stay_i %% 55L),
  age_interval = NA_character_,
  age_model_representation = "numeric",
  sex_class = ifelse(stay_i %% 2L == 0L, "female", "male"),
  cci_value = stay_i %% 8L,
  stringsAsFactors = FALSE)
topcoded <- stay_i %% 25L == 0L
stay_clinical$age_years_numeric[topcoded] <- NA_real_
stay_clinical$age_interval[topcoded] <- "90+"
stay_clinical$age_model_representation[topcoded] <- "source_interval"
write_tsv(stay_clinical, "stay_clinical.tsv")
units <- c(pao2 = "mmHg", paco2 = "mmHg", ph = "pH", lactate = "mmol/L", hb = "g/dL")
make_clinical <- function(pair, window) {
  n <- nrow(pair)
  i <- seq_len(n)
  patient_i <- match(pair$patient_key_internal, stays$patient_key_internal)
  x <- pair[ci1_pair_keys()]
  x$pair_key_internal <- sprintf("invented_rich_%s_%05d", window, i)
  x$sofa_total <- as.integer((patient_i + i) %% 18L)
  x$sofa_status <- "available"
  values <- list(
    pao2 = 57 + (patient_i %% 30L) + (i %% 9L) * 1.3,
    paco2 = 31 + (patient_i %% 21L) + (i %% 7L) * 0.9,
    ph = 7.18 + (patient_i %% 18L) * 0.016 + (i %% 5L) * 0.004,
    lactate = 0.8 + (patient_i %% 12L) * 0.28 + (i %% 4L) * 0.07,
    hb = 8.5 + (patient_i %% 15L) * 0.32 + (i %% 6L) * 0.08)
  for (name in names(units)) {
    missing <- (i + match(name, names(units))) %% 17L == 0L
    value <- values[[name]]
    value[missing] <- NA_real_
    x[[paste0(name, "_candidate_value")]] <- value
    x[[paste0(name, "_canonical_unit")]] <- ifelse(missing, "", units[[name]])
    x[[paste0(name, "_linkage_class")]] <- ifelse(missing, "", "same_specimen")
    x[[paste0(name, "_value_state")]] <- ifelse(missing, "unavailable",
                                                    "eligible_pending_candidate_build")
  }
  for (name in ci1_treatment_names()) {
    x[[paste0(name, "_state")]] <- ifelse(i %% 4L == 0L, "yes", "unknown")
    x[[paste0(name, "_scope")]] <- "permitted_by_gate_r08"
  }
  x[ci1_pair_columns()]
}
for (window in c(60L, 5L)) {
  matched <- r05_pair_cohort(loaded$spo2, loaded$sao2, "MIMIC", window)$final_pairs
  expected <- if (window == 60L) nrow(pairs$M60) else nrow(pairs$M5)
  if (nrow(matched) != expected) stop("RICH_SYNTHETIC_PAIR_COUNT_INVALID", call. = FALSE)
  write_tsv(make_clinical(matched, window), paste0("pair_clinical_", window, ".tsv"))
}
cat("RICH_SYNTHETIC_STANDARDIZED_INPUT_PASS\n")
