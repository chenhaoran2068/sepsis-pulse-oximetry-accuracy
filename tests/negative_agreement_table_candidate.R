#!/usr/bin/env Rscript
# Rejection tests operate on copied invented aggregate rows only.
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
write_tsv <- function(data, name) {
  path <- file.path(test_dir, paste0(name, ".tsv"))
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}
must_reject <- function(name, m, requested_window, marker) {
  mp <- write_tsv(m, paste0(name, "-manifest"))
  out <- file.path(test_dir, paste0(name, ".pdf"))
  rows <- file.path(test_dir, paste0(name, "-rows.tsv"))
  receipt <- file.path(test_dir, paste0(name, "-receipt.tsv"))
  log <- suppressWarnings(system2(rscript,
    c(shQuote(renderer), shQuote(mp), requested_window, shQuote(out), shQuote(rows), shQuote(receipt)),
    stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  if (status == 0L || !any(grepl(marker, log, fixed = TRUE)) ||
      any(file.exists(c(out, rows, receipt)))) {
    stop(paste("NEGATIVE_GUARD_FAILED", name, paste(log, collapse = " | ")), call. = FALSE)
  }
  cat("REJECTED", name, marker, "\n")
}
mutate_one <- function(name, change) {
  m <- manifest
  x <- read.delim(m$agreement_tsv[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
  x <- change(x)
  m$agreement_tsv[[1L]] <- write_tsv(x, paste0(name, "-aggregate"))
  m
}
bad <- manifest
bad$cohort[[2L]] <- bad$cohort[[1L]]
must_reject("duplicate-cohort", bad, window, "FIVE_COHORT_SET_INVALID")
must_reject("wrong-window", manifest, if (window == "60") "5" else "60",
            "COHORT_WINDOW_OR_ROW_MISMATCH")
must_reject("invalid-ci", mutate_one("invalid-ci", function(x) {
  x$mean_bias_ci_lower[[1L]] <- x$mean_bias[[1L]] + 0.01
  x
}), window, "AGREEMENT_CI_INVALID")
must_reject("reversed-loa", mutate_one("reversed-loa", function(x) {
  x$lower_loa[[1L]] <- x$mean_bias[[1L]] + 0.01
  x$lower_loa_ci_lower[[1L]] <- x$lower_loa[[1L]] - 0.2
  x$lower_loa_ci_upper[[1L]] <- x$lower_loa[[1L]] + 0.2
  x
}), window, "LOA_ORDER_INVALID")
cat("AGREEMENT_TABLE_NEGATIVE_QA_PASS checks=4\n")
