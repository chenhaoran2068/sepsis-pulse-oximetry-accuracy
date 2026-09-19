#!/usr/bin/env Rscript
# Malformed entirely invented aggregate inputs must fail before display creation.
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
  p <- file.path(test_dir, paste0(name, ".tsv"))
  write.table(data, p, row.names = FALSE, quote = FALSE, sep = "\t")
  p
}
must_reject <- function(name, manifest, marker) {
  mp <- write_tsv(manifest, paste0(name, "-manifest"))
  pdf <- file.path(test_dir, paste0(name, ".pdf"))
  rows <- file.path(test_dir, paste0(name, "-rows.tsv"))
  receipt <- file.path(test_dir, paste0(name, "-receipt.tsv"))
  out <- suppressWarnings(system2(rscript,
    c(shQuote(renderer), shQuote(mp), shQuote(pdf), shQuote(rows), shQuote(receipt)),
    stdout = TRUE, stderr = TRUE))
  code <- attr(out, "status")
  if (is.null(code)) code <- 0L
  if (code == 0L || !any(grepl(marker, out, fixed = TRUE)) ||
      any(file.exists(c(pdf, rows, receipt)))) {
    stop(paste("NEGATIVE_GUARD_FAILED", name, paste(out, collapse = " | ")), call. = FALSE)
  }
  cat("REJECTED", name, marker, "\n")
}
mutate_one <- function(name, change) {
  copy <- m
  x <- read.delim(copy$low_sao2_tsv[[1L]], check.names = FALSE, stringsAsFactors = FALSE)
  x <- change(x)
  copy$low_sao2_tsv[[1L]] <- write_tsv(x, paste0(name, "-aggregate"))
  copy
}
bad <- m
bad$window_minutes[[2L]] <- 60L
must_reject("duplicate-window", bad, "TEN_COHORT_WINDOW_SET_INVALID")
must_reject("non-nested-count", mutate_one("non-nested-count", function(x) {
  x$pair_numerator_n[x$threshold_spo2_ge == 92] <- x$pair_numerator_n[x$threshold_spo2_ge == 88] + 1L
  x$pair_conditional_proportion[x$threshold_spo2_ge == 92] <-
    x$pair_numerator_n[x$threshold_spo2_ge == 92] /
    x$pair_denominator_sao2_lt88_n[x$threshold_spo2_ge == 92]
  x$pair_ci_upper[x$threshold_spo2_ge == 92] <- 0.99
  x$pair_ci_lower[x$threshold_spo2_ge == 92] <- 0.1
  x
}), "THRESHOLD_NESTING_INVALID")
must_reject("bad-proportion", mutate_one("bad-proportion", function(x) {
  x$pair_conditional_proportion[[1L]] <- 0.1
  x
}), "THRESHOLD_ESTIMATE_OR_CI_INVALID")
must_reject("bad-denominator", mutate_one("bad-denominator", function(x) {
  x$pair_denominator_sao2_lt88_n[[1L]] <- x$accepted_pair_n[[1L]] + 1L
  x
}), "THRESHOLD_COUNTS_OR_STATUS_INVALID")
sparse <- mutate_one("valid-sparse", function(x) {
  x$cell_status[[1L]] <- "SPARSE"
  x$pair_conditional_proportion[[1L]] <- NA_real_
  x$pair_ci_lower[[1L]] <- NA_real_
  x$pair_ci_upper[[1L]] <- NA_real_
  x
})
sp <- write_tsv(sparse, "valid-sparse-manifest")
pdf <- file.path(test_dir, "valid-sparse.pdf")
rows <- file.path(test_dir, "valid-sparse-rows.tsv")
receipt <- file.path(test_dir, "valid-sparse-receipt.tsv")
out <- suppressWarnings(system2(rscript,
  c(shQuote(renderer), shQuote(sp), shQuote(pdf), shQuote(rows), shQuote(receipt)),
  stdout = TRUE, stderr = TRUE))
code <- attr(out, "status")
if ((!is.null(code) && code != 0L) || !all(file.exists(c(pdf, rows, receipt)))) {
  stop(paste("VALID_SPARSE_REJECTED", paste(out, collapse = " | ")), call. = FALSE)
}
display <- read.delim(rows, check.names = FALSE, stringsAsFactors = FALSE)
if (grepl("%", display$occult_n_N_pct[[1L]], fixed = TRUE) ||
    grepl("(", display$occult_n_N_pct[[1L]], fixed = TRUE) ||
    display$occult_ci_pct[[1L]] != "-") {
  stop("VALID_SPARSE_DISPLAY_INVALID", call. = FALSE)
}
cat("THRESHOLD_TABLE_GUARD_QA_PASS checks=5\n")
