#!/usr/bin/env Rscript
# Independent rejection checks, using copies of entirely invented fixture rows.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: VALID_MANIFEST.tsv WINDOW RENDERER.R NEW_TEST_DIRECTORY", call. = FALSE)
manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
window <- args[[2L]]
renderer <- normalizePath(args[[3L]], mustWork = TRUE)
test_dir <- args[[4L]]
if (file.exists(test_dir)) stop("TEST_DIRECTORY_EXISTS", call. = FALSE)
dir.create(test_dir, recursive = TRUE)
manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
write_tsv <- function(data, stem) {
  path <- file.path(test_dir, paste0(stem, ".tsv"))
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}
must_reject <- function(stem, m, requested_window, marker) {
  mp <- write_tsv(m, paste0(stem, "-manifest"))
  pdf <- file.path(test_dir, paste0(stem, ".pdf"))
  receipt <- file.path(test_dir, paste0(stem, "-receipt.tsv"))
  output <- suppressWarnings(system2(rscript,
    c(shQuote(renderer), shQuote(mp), requested_window, shQuote(pdf), shQuote(receipt)),
    stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (status == 0L || !any(grepl(marker, output, fixed = TRUE)) ||
      file.exists(pdf) || file.exists(receipt)) {
    stop(paste("NEGATIVE_GUARD_FAILED", stem, paste(output, collapse = " | ")), call. = FALSE)
  }
  cat("REJECTED", stem, marker, "\n")
}
mutate_one <- function(stem, change) {
  m <- manifest
  x <- read.delim(m$day_tsv[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
  x <- change(x)
  m$day_tsv[[1L]] <- write_tsv(x, paste0(stem, "-day"))
  m
}
bad <- manifest
bad$cohort[[2L]] <- bad$cohort[[1L]]
must_reject("duplicate-cohort", bad, window, "FIVE_COHORT_SET_INVALID")
other_window <- if (window == "60") "5" else "60"
must_reject("wrong-window", manifest, other_window, "COHORT_WINDOW_OR_CELL_COUNT_MISMATCH")
must_reject("missing-grid-cell", mutate_one("missing-grid-cell", function(x) {
  x$spo2_stratum[[1L]] <- x$spo2_stratum[[2L]]
  x
}), window, "DAY_CELL_GRID_INVALID")
must_reject("impossible-count", mutate_one("impossible-count", function(x) {
  x$numerator_n[[1L]] <- x$denominator_n[[1L]] + 1L
  x
}), window, "DAY_COUNT_OR_STATUS_INVALID")
must_reject("sparse-marked-displayable", mutate_one("sparse-marked-displayable", function(x) {
  x$patient_n[[1L]] <- 9L
  x
}), window, "DISPLAYABLE_CELL_IS_SPARSE")
must_reject("bad-ci", mutate_one("bad-ci", function(x) {
  x$ci_lower[[1L]] <- x$conditional_proportion[[1L]] + 0.01
  x
}), window, "DAY_ESTIMATE_OR_CI_INVALID")
cat("ICU_DAY_NEGATIVE_QA_PASS checks=6\n")
