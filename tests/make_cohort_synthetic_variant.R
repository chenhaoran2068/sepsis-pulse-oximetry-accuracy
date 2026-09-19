args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: RICH_MIMIC_STANDARDIZED_DIR COHORT NEW_OUTPUT_DIR", call. = FALSE)
input <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cohort <- args[[2L]]
output <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
if (!cohort %in% c("Amsterdam", "eICU", "SICDB", "Lianyungang"))
  stop("COHORT_VARIANT_INVALID", call. = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
required <- c("stays.tsv", "spo2_events.tsv", "sao2_events.tsv", "stay_clinical.tsv",
              "pair_clinical_60.tsv", "pair_clinical_5.tsv")
if (!all(file.exists(file.path(input, required)))) stop("INPUT_FILES_MISSING", call. = FALSE)
read_tsv <- function(name) read.delim(file.path(input, name), check.names = FALSE,
                                     stringsAsFactors = FALSE, na.strings = "")
write_tsv <- function(x, name) utils::write.table(x, file.path(output, name), sep = "\t",
                                                   row.names = FALSE, quote = FALSE, na = "")
dir.create(output, recursive = TRUE)
source_stays <- read_tsv("stays.tsv")
patient_for_stay <- if (cohort == "eICU")
  sprintf("invented_eicu_patient_%04d", ceiling(seq_len(nrow(source_stays)) / 2)) else
  source_stays$patient_key_internal
translate_patient <- function(x) {
  if (cohort != "eICU") return(x)
  index <- match(x$stay_key_internal, source_stays$stay_key_internal)
  if (anyNA(index)) stop("STAY_ID_NOT_FOUND", call. = FALSE)
  x$patient_key_internal <- patient_for_stay[index]
  x
}
for (name in c("stays.tsv", "spo2_events.tsv", "sao2_events.tsv"))
  write_tsv(translate_patient(read_tsv(name)), name)
stays <- read_tsv("stay_clinical.tsv")
if (cohort == "eICU") {
  first <- rep(seq(1L, nrow(stays), by = 2L), each = 2L)[seq_len(nrow(stays))]
  for (field in c("age_years_numeric", "age_interval", "age_model_representation",
                  "sex_class", "cci_value"))
    stays[[field]] <- stays[[field]][first]
  # The source fixture alternates sex by stay. Selecting every first stay would
  # otherwise create an all-male invented eICU sample with an unfittable factor.
  stays$sex_class <- ifelse(ceiling(seq_len(nrow(stays)) / 2) %% 2L == 0L,
                            "female", "male")
  stays <- translate_patient(stays)
}
age <- stays$age_years_numeric
top <- is.na(age) & stays$age_interval == "90+"
if (!any(top) || any(is.na(age) & !top)) stop("SOURCE_AGE_FIXTURE_INVALID", call. = FALSE)
if (cohort == "Amsterdam") {
  starts <- c(18, 40, 50, 60, 70, 80)
  labels <- c("18-39", "40-49", "50-59", "60-69", "70-79", "80+")
  age[top] <- 90
  stays$age_interval <- labels[findInterval(age, starts)]
  stays$age_years_numeric <- NA_real_
  stays$age_model_representation <- "source_interval"
} else if (cohort == "SICDB") {
  age[top] <- 90
  rounded <- 5 * floor((age + 2.5) / 5)
  rounded <- pmax(20, pmin(85, rounded))
  # The formal factor contract has 20- and 25-year bins. Populate them in this
  # invented fixture so empty design columns do not masquerade as a code defect.
  ordinary <- which(age < 90)
  rounded[ordinary[seq_len(10L)]] <- 20
  rounded[ordinary[11:20]] <- 25
  stays$age_interval <- ifelse(age >= 90, "90+_rounded_or_topcoded",
                              paste0(rounded, "_rounded_5y_bin"))
  stays$age_years_numeric <- NA_real_
  stays$age_model_representation <- "source_interval"
} else if (cohort == "eICU") {
  stays$age_interval[top] <- "90+"
} else {
  stays$age_years_numeric[top] <- 90
  stays$age_interval <- NA_character_
  stays$age_model_representation <- "numeric"
}
write_tsv(stays, "stay_clinical.tsv")
for (window in c(60L, 5L)) {
  pairs <- read_tsv(paste0("pair_clinical_", window, ".tsv"))
  pairs <- translate_patient(pairs)
  if (cohort %in% c("eICU", "SICDB", "Lianyungang")) {
    # Keep some PaCO2 values missing to exercise MI, but avoid making this
    # invented formal-fit fixture hinge on Monte Carlo noise at a 0.10 gate.
    absent <- which(is.na(pairs$paco2_candidate_value))
    fill <- absent[seq_along(absent) %% 5L != 0L]
    pairs$paco2_candidate_value[fill] <- 38 + (fill %% 15L) / 2
    pairs$paco2_canonical_unit[fill] <- "mmHg"
    pairs$paco2_linkage_class[fill] <- "same_specimen"
    pairs$paco2_value_state[fill] <- "eligible_pending_candidate_build"
  }
  if (cohort %in% c("eICU", "SICDB", "Lianyungang")) {
    pairs$current_invasive_mechanical_ventilation_state <- NA_character_
    pairs$current_invasive_mechanical_ventilation_scope <- "excluded_by_gate_r08"
  }
  if (cohort == "Lianyungang") {
    pairs$concurrent_vasoactive_medication_use_state <- NA_character_
    pairs$concurrent_vasoactive_medication_use_scope <- "excluded_by_gate_r08"
  }
  write_tsv(pairs, paste0("pair_clinical_", window, ".tsv"))
}
cat("COHORT_SYNTHETIC_VARIANT_PASS cohort=", cohort, "\n")
