#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Usage: ANALYSIS_MANIFEST.tsv ROWS.tsv RECEIPT.tsv",call.=FALSE)
m <- fread(args[[1L]])
rows <- fread(args[[2L]])
receipt <- fread(args[[3L]])
if (nrow(m)!=5L || nrow(rows)!=100L || nrow(receipt)!=1L ||
    receipt$row_n!=100L || anyDuplicated(rows[,.(cohort,window_minutes,characteristic)]))
  stop("PAIR_CHARACTERISTICS_QA_SHAPE_FAIL",call.=FALSE)
display_names <- c(MIMIC="MIMIC",Amsterdam="AmsterdamUMCdb",eICU="eICU",
                   SICDB="SICdb",Lianyungang="Lianyungang")
checks <- 0L
for (i in seq_len(nrow(m))) for (w in c(60L,5L)) {
  p <- as.data.table(readRDS(file.path(m$analysis_dir[[i]],
                                    paste0("matched_pairs_",w,".rds"))))
  r <- rows[cohort==display_names[[m$cohort[[i]]]] & window_minutes==w]
  n <- nrow(p)
  get <- function(label) r[characteristic==label]$display_value
  if (nrow(r)!=10L || as.integer(gsub(",","",get("Eligible pairs, n")))!=n)
    stop("PAIR_CHARACTERISTICS_QA_COUNT_FAIL",call.=FALSE)
  counts <- c(sum(p$pair_direction=="spo2_before_sao2"),
              sum(p$pair_direction=="same_time"),
              sum(p$pair_direction=="spo2_after_sao2"))
  labels <- c("SpO2 before SaO2, n (%)","SpO2 simultaneous with SaO2, n (%)",
              "SpO2 after SaO2, n (%)")
  for (j in seq_along(labels)) {
    displayed <- as.integer(gsub(",","",sub(" .*","",get(labels[[j]]))))
    if (displayed!=counts[[j]]) stop("PAIR_CHARACTERISTICS_QA_DIRECTION_FAIL",call.=FALSE)
    checks <- checks+1L
  }
  values <- c(mean(p$spo2_saturation_percent),mean(p$sao2_saturation_percent),
              mean((p$spo2_saturation_percent+p$sao2_saturation_percent)/2))
  measures <- c("SpO2 in eligible pairs, %, mean (SD)",
                "SaO2 in eligible pairs, %, mean (SD)",
                "Paired mean saturation, %, mean (SD)")
  for (j in seq_along(measures)) {
    displayed <- as.numeric(sub(" .*","",get(measures[[j]])))
    if (abs(displayed-values[[j]])>.051)
      stop("PAIR_CHARACTERISTICS_QA_MEAN_FAIL",call.=FALSE)
    checks <- checks+1L
  }
  checks <- checks+2L
}
cat("PAIR_CHARACTERISTICS_INDEPENDENT_QA_PASS checks=",checks,"\n",sep="")
