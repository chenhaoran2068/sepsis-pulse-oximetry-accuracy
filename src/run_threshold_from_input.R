args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L || !args[[5L]] %in% c("paper", "smoke"))
  stop("Usage: COHORT WINDOW_MINUTES ANALYSIS_INPUT_DIR NEW_OUTPUT_DIR paper|smoke", call. = FALSE)
cohort <- args[[1L]]
window <- as.integer(args[[2L]])
input <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[4L]], winslash = "/", mustWork = FALSE)
mode <- args[[5L]]
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
if (!cohort %in% cohorts || !window %in% c(60L, 5L)) stop("THRESHOLD_SCOPE_INVALID", call. = FALSE)
if (file.exists(output) || dir.exists(output)) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)
for (name in c("stage404_core.R", "stage404b_core.R", "stage404c_core.R"))
  source(file.path(root, "code", "cores", name))
suppressPackageStartupMessages({library(data.table); library(digest)})
pair_path <- file.path(input, paste0("matched_pairs_", window, ".rds"))
analytic_path <- file.path(input, paste0("analytic_pairs_", window, ".rds"))
day_path <- file.path(input, paste0("day1_7_pairs_", window, ".rds"))
receipt <- read.delim(file.path(input, "receipt.tsv"), check.names = FALSE)
hash_ok <- function(path) {
  row <- receipt[receipt$file == basename(path), ]
  file.exists(path) && nrow(row) == 1L &&
    tolower(row$sha256[[1L]]) == tolower(digest::digest(file = path, algo = "sha256"))
}
if (!all(vapply(c(pair_path, analytic_path, day_path), hash_ok, logical(1))))
  stop("THRESHOLD_INPUT_HASH_INVALID", call. = FALSE)
full <- data.table::as.data.table(readRDS(pair_path))
analytic <- data.table::as.data.table(readRDS(analytic_path))
if (nrow(full) != nrow(analytic) ||
    !identical(as.character(full$stay_key_internal), as.character(analytic$stay_key_internal)) ||
    !identical(as.character(full$patient_key_internal), as.character(analytic$patient_key_internal)))
  stop("THRESHOLD_PAIR_IDENTITY_INVALID", call. = FALSE)
full[, pair_key_internal := analytic$pair_key_internal]
window_label <- if (window == 60L) "M60" else "M5"
s404_validate_full_pairs(full, cohort, window_label)
q21 <- data.table::as.data.table(readRDS(day_path))
s404_validate_q21_pairs(q21, cohort, window_label)
q21[, display_stratum := s404_source_to_display_stratum(spo2_display_stratum)]
if (anyNA(q21$display_stratum)) stop("THRESHOLD_DAY_STRATUM_MAP_INVALID", call. = FALSE)
reps <- if (mode == "paper") 2000L else 30L
cohort_index <- match(cohort, cohorts)
window_index <- if (window == 60L) 1L else 2L
dir.create(output, recursive = TRUE)
writeLines("THRESHOLD_RUNNING", file.path(output, "status.txt"))
write_tsv <- function(x, name) data.table::fwrite(x, file.path(output, name), sep = "\t")
tryCatch({
  forward <- rbindlist(lapply(seq_along(c(88, 92)), function(i) {
    threshold <- c(88, 92)[[i]]
    s404_overall_pattern(full, threshold, reps,
      404000L + window_index * 10000L + cohort_index * 100L + i)
  }))
  reverse <- rbindlist(lapply(seq_along(c(88, 92)), function(i) {
    threshold <- c(88, 92)[[i]]
    s404b_reverse_pattern(full, threshold, reps,
      404500L + window_index * 10000L + cohort_index * 100L + i)
  }))
  patients <- rbindlist(lapply(c(88, 92), function(threshold)
    s404c_pattern_summary(copy(full), threshold)$patient_summary))
  concentration <- rbindlist(lapply(c(88, 92), function(threshold)
    s404c_pattern_summary(copy(full), threshold)$concentration))
  cells <- list()
  position <- 0L
  for (day in 1:7) for (i in seq_along(s404_strata())) {
    position <- position + 1L
    stratum <- s404_strata()[[i]]
    subset <- q21[relative_icu_day == day & display_stratum == stratum]
    seed <- 409000L + window_index * 10000L + cohort_index * 1000L + day * 10L + i
    row <- s404_day_stratum_cell(subset, reps, seed)
    row[, `:=`(cohort = cohort, window_minutes = window, relative_icu_day = day,
               spo2_stratum = stratum,
               reportable_percentage = ifelse(cell_status == "DISPLAYABLE",
                                              100 * conditional_proportion, NA_real_))]
    cells[[position]] <- row
  }
  write_tsv(forward, "threshold_from_low_sao2.tsv")
  write_tsv(reverse, "threshold_from_displayed_spo2.tsv")
  write_tsv(patients, "threshold_patient_summary.tsv")
  write_tsv(concentration, "threshold_event_concentration.tsv")
  write_tsv(rbindlist(cells), "threshold_day1_7_strata.tsv")
  write_tsv(data.table(cohort = cohort, window_minutes = window, mode = mode,
                       full_pair_n = nrow(full), day1_7_pair_n = nrow(q21),
                       bootstrap_replicates = reps,
                       full_sha256 = digest::digest(file = pair_path, algo = "sha256"),
                       day_sha256 = digest::digest(file = day_path, algo = "sha256")),
            "run_receipt.tsv")
  writeLines(if (mode == "paper") "THRESHOLD_PAPER_SETTINGS_COMPLETE_CANDIDATE_ONLY" else
               "THRESHOLD_SMOKE_COMPLETE_NOT_PAPER_SETTINGS", file.path(output, "status.txt"))
  cat("THRESHOLD_RUN_PASS mode=", mode, "\n")
}, error = function(e) {
  writeLines(paste("THRESHOLD_FAIL", conditionMessage(e)), file.path(output, "status.txt"))
  stop(e)
})
