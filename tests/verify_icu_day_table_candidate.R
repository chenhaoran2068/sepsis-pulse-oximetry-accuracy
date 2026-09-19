#!/usr/bin/env Rscript
# Independent 140-cell reproduction check of the two-page ICU-day table.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: MANIFEST.tsv WINDOW OUTPUT.pdf ROWS.tsv RECEIPT.tsv", call. = FALSE)
m <- fread(args[[1L]])
window <- as.integer(args[[2L]])
rows <- fread(args[[4L]], colClasses = "character")
receipt <- fread(args[[5L]])
checks <- character()
check <- function(name, ok) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(paste("QA_FAIL", name), call. = FALSE)
  checks <<- c(checks, name)
}
check("pdf_magic", identical(readBin(args[[3L]], "raw", n = 4L), charToRaw("%PDF")))
check("rows", nrow(rows) == 35L && ncol(rows) == 6L)
check("receipt", nrow(receipt) == 1L && receipt$window_minutes[[1L]] == window &&
      receipt$input_cell_n[[1L]] == 140L && receipt$display_row_n[[1L]] == 35L &&
      receipt$page_n[[1L]] == 2L)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
labels <- c("MIMIC", "AmsterdamUMCdb", "eICU", "SICdb", "Lianyungang")
strata <- c("70%-<88%", "88%-<92%", "92%-<=96%", ">96%-100%")
fields <- c("spo2_70_to_lt88", "spo2_88_to_lt92", "spo2_92_to_le96", "spo2_gt96_to_100")
fmt_n <- function(z) format(z, big.mark = ",", scientific = FALSE, trim = TRUE)
sparse_total <- 0L
for (i in seq_along(cohorts)) {
  input <- fread(m[cohort == cohorts[[i]], day_tsv])
  check(paste0(cohorts[[i]], "_input_shape"), nrow(input) == 28L &&
        all(input$window_minutes == window))
  sparse_total <- sparse_total + sum(input$cell_status != "DISPLAYABLE")
  for (d in 1:7) {
    r <- rows[(i - 1L) * 7L + d]
    check(paste0("row_", i, "_", d),
          identical(r$cohort[[1L]], if (d == 1L) labels[[i]] else "") &&
          identical(r$icu_day[[1L]], as.character(d)))
    for (s in seq_along(strata)) {
      x <- input[relative_icu_day == d & spo2_stratum == strata[[s]]]
      check(paste0("input_cell_", i, "_", d, "_", s), nrow(x) == 1L)
      count <- paste0(fmt_n(x$numerator_n), "/", fmt_n(x$denominator_n))
      expected <- if (x$cell_status == "DISPLAYABLE") {
        sprintf("%s (%.1f) [%.1f to %.1f]", count, 100 * x$conditional_proportion,
                100 * x$ci_lower, 100 * x$ci_upper)
      } else count
      check(paste0("cell_", i, "_", d, "_", s), identical(r[[fields[[s]]]][[1L]], expected))
    }
  }
}
check("sparse_receipt", receipt$sparse_n[[1L]] == sparse_total)
cat("ICU_DAY_TABLE_INDEPENDENT_QA_PASS window=", window,
    " checks=", length(checks), "\n", sep = "")
