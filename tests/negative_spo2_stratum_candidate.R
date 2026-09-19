#!/usr/bin/env Rscript
# Guard tests use only copied synthetic aggregate inputs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: VALID_MANIFEST.tsv RENDERER.R NEW_TEST_DIRECTORY", call. = FALSE)
manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
renderer <- normalizePath(args[[2L]], mustWork = TRUE)
test_dir <- args[[3L]]
if (file.exists(test_dir)) stop("TEST_DIRECTORY_EXISTS", call. = FALSE)
dir.create(test_dir, recursive = TRUE)
m <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
write_tsv <- function(data, name) {
  path <- file.path(test_dir, paste0(name, ".tsv"))
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}
must_reject <- function(name, manifest, marker) {
  mp <- write_tsv(manifest, paste0(name, "-manifest"))
  pdf <- file.path(test_dir, paste0(name, ".pdf"))
  receipt <- file.path(test_dir, paste0(name, "-receipt.tsv"))
  output <- suppressWarnings(system2(rscript,
    c(shQuote(renderer), shQuote(mp), shQuote(pdf), shQuote(receipt)),
    stdout = TRUE, stderr = TRUE))
  code <- attr(output, "status")
  if (is.null(code)) code <- 0L
  if (code == 0L || !any(grepl(marker, output, fixed = TRUE)) ||
      file.exists(pdf) || file.exists(receipt)) {
    stop(paste("NEGATIVE_GUARD_FAILED", name, paste(output, collapse = " | ")), call. = FALSE)
  }
  cat("REJECTED", name, marker, "\n")
}
bad <- m
bad$cohort[[2L]] <- bad$cohort[[1L]]
must_reject("duplicate-cohort", bad, "FIVE_COHORT_SET_INVALID")
# The renderer derives the stratum file from the aggregate's parent directory.
# Create an isolated parent for the manipulated pair of files.
mutate_one <- function(name, change) {
  parent <- file.path(test_dir, name)
  dir.create(parent)
  original <- m$agreement_tsv[[1L]]
  file.copy(original, file.path(parent, "agreement_overall.tsv"))
  x <- read.delim(file.path(dirname(original), "agreement_by_spo2_stratum.tsv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
  write.table(change(x), file.path(parent, "agreement_by_spo2_stratum.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  copy <- m
  copy$agreement_tsv[[1L]] <- file.path(parent, "agreement_overall.tsv")
  copy
}
must_reject("wrong-window", mutate_one("wrong-window", function(x) {
  x$window_minutes[[1L]] <- 5L
  x
}), "STRATUM_AGREEMENT_GRID_INVALID")
must_reject("invalid-loa", mutate_one("invalid-loa", function(x) {
  x$lower_loa[[1L]] <- x$mean_bias[[1L]] + 1
  x
}), "STRATUM_AGREEMENT_VALUES_INVALID")
cat("STRATUM_AGREEMENT_NEGATIVE_QA_PASS checks=3\n")
