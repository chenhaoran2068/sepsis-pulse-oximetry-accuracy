#!/usr/bin/env Rscript
# Guard tests operate only on a new synthetic test directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: VALID_MANIFEST.tsv RENDERER.R NEW_TEST_DIRECTORY", call. = FALSE)
valid_manifest <- normalizePath(args[[1L]], mustWork = TRUE)
renderer <- normalizePath(args[[2L]], mustWork = TRUE)
test_dir <- args[[3L]]
if (file.exists(test_dir)) stop("TEST_DIRECTORY_EXISTS", call. = FALSE)
dir.create(test_dir, recursive = TRUE)
manifest <- read.delim(valid_manifest, check.names = FALSE, stringsAsFactors = FALSE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")

write_manifest <- function(data, name) {
  path <- file.path(test_dir, name)
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}
must_reject <- function(test_name, manifest_path, window, expected_error) {
  pdf <- file.path(test_dir, paste0(test_name, ".pdf"))
  receipt <- file.path(test_dir, paste0(test_name, ".tsv"))
  output <- suppressWarnings(system2(rscript,
    c(shQuote(renderer), shQuote(manifest_path), window, shQuote(pdf), shQuote(receipt)),
    stdout = TRUE, stderr = TRUE))
  code <- attr(output, "status")
  if (is.null(code)) code <- 0L
  if (code == 0L || !any(grepl(expected_error, output, fixed = TRUE)) ||
      file.exists(pdf) || file.exists(receipt)) {
    stop(paste("NEGATIVE_GUARD_FAILED", test_name, paste(output, collapse = " | ")), call. = FALSE)
  }
  cat("REJECTED", test_name, expected_error, "\n")
}

bad_cohort <- manifest
bad_cohort$cohort[[2L]] <- "MIMIC"
must_reject("duplicate-cohort", write_manifest(bad_cohort, "duplicate-cohort-manifest.tsv"),
            "60", "FIVE_COHORT_SET_INVALID")
must_reject("wrong-window", valid_manifest, "5", "COHORT_WINDOW_OR_COUNT_MISMATCH")
bad_summary <- manifest
summary_copy <- read.delim(manifest$agreement_tsv[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
summary_copy$pair_n[[1L]] <- summary_copy$pair_n[[1L]] + 1L
summary_path <- file.path(test_dir, "wrong-count-summary.tsv")
write.table(summary_copy, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
bad_summary$agreement_tsv[[1L]] <- summary_path
must_reject("wrong-count", write_manifest(bad_summary, "wrong-count-manifest.tsv"),
            "60", "COHORT_WINDOW_OR_COUNT_MISMATCH")
cat("BLAND_ALTMAN_NEGATIVE_QA_PASS checks=3\n")
