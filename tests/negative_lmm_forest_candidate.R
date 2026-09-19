#!/usr/bin/env Rscript
# Reject malformed pooled coefficients before figure generation.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: VALID_MANIFEST.tsv WINDOW RENDERER.R NEW_TEST_DIRECTORY", call. = FALSE)
manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
window <- args[[2L]]
renderer <- normalizePath(args[[3L]], mustWork = TRUE)
test_dir <- args[[4L]]
if (file.exists(test_dir)) stop("TEST_DIRECTORY_EXISTS", call. = FALSE)
dir.create(test_dir, recursive = TRUE)
m <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
write_tsv <- function(x, name) {
  p <- file.path(test_dir, paste0(name, ".tsv"))
  write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE)
  p
}
must_reject <- function(name, manifest, marker) {
  mp <- write_tsv(manifest, paste0(name, "-manifest"))
  pdf <- file.path(test_dir, paste0(name, ".pdf"))
  chart <- file.path(test_dir, paste0(name, "-chart.tsv"))
  receipt <- file.path(test_dir, paste0(name, "-receipt.tsv"))
  log <- suppressWarnings(system2(rscript,
    c(shQuote(renderer), shQuote(mp), window, shQuote(pdf), shQuote(chart), shQuote(receipt)),
    stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  if (status == 0L || !any(grepl(marker, log, fixed = TRUE)) ||
      any(file.exists(c(pdf, chart, receipt)))) {
    stop(paste("NEGATIVE_GUARD_FAILED", name, paste(log, collapse = " | ")), call. = FALSE)
  }
  cat("REJECTED", name, marker, "\n")
}
mutate_one <- function(name, change) {
  copy <- m
  x <- read.delim(copy$lmm_tsv[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
  copy$lmm_tsv[[1L]] <- write_tsv(change(x), paste0(name, "-pooled"))
  copy
}
bad <- m
bad$cohort[[2L]] <- bad$cohort[[1L]]
must_reject("duplicate-cohort", bad, "FIVE_COHORT_SET_INVALID")
must_reject("duplicate-term", mutate_one("duplicate-term", function(x) {
  rbind(x, x[1L, , drop = FALSE])
}), "LMM_DUPLICATE_TERM")
must_reject("invalid-multiplier", mutate_one("invalid-multiplier", function(x) {
  x$reporting_multiplier[[1L]] <- 0
  x
}), "LMM_INPUT_VALUE_INVALID")
must_reject("invalid-ci", mutate_one("invalid-ci", function(x) {
  x$ci_lower[[2L]] <- x$estimate[[2L]] + 1
  x
}), "LMM_CI_INVALID")
must_reject("missing-paired-mean", mutate_one("missing-paired-mean", function(x) {
  x[x$term != "paired_mean_saturation_c90", , drop = FALSE]
}), "PAIRED_MEAN_MISSING")
cat("LMM_FOREST_NEGATIVE_QA_PASS checks=5\n")
