#!/usr/bin/env Rscript
# Technical renderer for Supplementary Table 10 from both independently paired sets.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: ANALYSIS_MANIFEST.tsv NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv", call. = FALSE)
cohorts <- c("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
pretty <- c(MIMIC="MIMIC", Amsterdam="AmsterdamUMCdb", eICU="eICU",
            SICDB="SICdb", Lianyungang="Lianyungang")
m <- fread(args[[1L]])
if (!identical(names(m), c("cohort", "analysis_dir")) || nrow(m) != 5L ||
    !setequal(m$cohort, cohorts) || anyDuplicated(m$cohort))
  stop("PAIR_CHARACTERISTICS_MANIFEST_INVALID", call. = FALSE)
out <- args[2:4]
if (any(file.exists(out))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
if (any(!dir.exists(dirname(out)))) stop("OUTPUT_PARENT_MISSING", call. = FALSE)
fmt_n <- function(x) format(x, big.mark=",", scientific=FALSE, trim=TRUE)
fmt_mean <- function(x) sprintf("%.1f (%.1f)", mean(x), sd(x))
fmt_median <- function(x) {
  q <- quantile(x, c(.25,.5,.75), names=FALSE, type=7)
  sprintf("%.1f [%.1f, %.1f]", q[[2L]], q[[1L]], q[[3L]])
}
fmt_count <- function(n, denominator) {
  pct <- 100 * n / denominator
  paste0(fmt_n(n), " (", if (n > 0L && pct < .05) "<0.1" else sprintf("%.1f", pct), ")")
}
result <- list()
for (w in c(60L,5L)) for (cohort_value in cohorts) {
  dir <- m[cohort == cohort_value]$analysis_dir
  path <- file.path(dir, paste0("matched_pairs_",w,".rds"))
  analytic_path <- file.path(dir, paste0("analytic_pairs_",w,".rds"))
  if (length(dir) != 1L || !file.exists(path) || !file.exists(analytic_path))
    stop("PAIR_CHARACTERISTICS_SOURCE_MISSING", call. = FALSE)
  p <- as.data.table(readRDS(path))
  analytic <- as.data.table(readRDS(analytic_path))
  need <- c("patient_key_internal","stay_key_internal","spo2_saturation_percent",
            "sao2_saturation_percent","absolute_lag_minutes","pair_direction")
  if (!all(need %in% names(p)) || nrow(p)<1L || nrow(p)!=nrow(analytic) ||
      anyNA(p[, ..need]) ||
      any(!is.finite(as.matrix(p[, .(spo2_saturation_percent,sao2_saturation_percent,
                                     absolute_lag_minutes)]))) ||
      any(p$absolute_lag_minutes < 0 | p$absolute_lag_minutes > w) ||
      any(!p$pair_direction %in%
          c("spo2_before_sao2","same_time","spo2_after_sao2")))
    stop("PAIR_CHARACTERISTICS_SOURCE_INVALID", call. = FALSE)
  pair_uid <- paste(p$patient_key_internal,p$stay_key_internal,
                    p$spo2_stable_event_ordinal,p$sao2_stable_event_ordinal,sep="|")
  analytic_uid <- paste(analytic$patient_key_internal,analytic$stay_key_internal,
                        analytic$spo2_stable_event_ordinal,analytic$sao2_stable_event_ordinal,sep="|")
  if (anyDuplicated(pair_uid) || !setequal(pair_uid,analytic_uid))
    stop("PAIR_CHARACTERISTICS_ACCEPTED_SET_MISMATCH", call. = FALSE)
  n <- nrow(p)
  patient_counts <- p[, .N, by=patient_key_internal]$N
  stay_counts <- p[, .N, by=.(patient_key_internal,stay_key_internal)]$N
  directions <- c(before="spo2_before_sao2",
                  simultaneous="same_time",after="spo2_after_sao2")
  vals <- c(
    "Eligible pairs, n"=fmt_n(n),
    "Pairs per patient, median [IQR]"=fmt_median(patient_counts),
    "Pairs per ICU stay, median [IQR]"=fmt_median(stay_counts),
    "SpO2 in eligible pairs, %, mean (SD)"=fmt_mean(p$spo2_saturation_percent),
    "SaO2 in eligible pairs, %, mean (SD)"=fmt_mean(p$sao2_saturation_percent),
    "Paired mean saturation, %, mean (SD)"=
      fmt_mean((p$spo2_saturation_percent+p$sao2_saturation_percent)/2),
    "Absolute time difference, min, median [IQR]"=fmt_median(p$absolute_lag_minutes),
    "SpO2 before SaO2, n (%)"=fmt_count(sum(p$pair_direction==directions[["before"]]),n),
    "SpO2 simultaneous with SaO2, n (%)"=fmt_count(sum(p$pair_direction==directions[["simultaneous"]]),n),
    "SpO2 after SaO2, n (%)"=fmt_count(sum(p$pair_direction==directions[["after"]]),n))
  result[[paste(cohort_value,w)]] <- data.table(cohort=unname(pretty[[cohort_value]]),
                                                window_minutes=w, characteristic=names(vals),
                                                display_value=unname(vals))
}
rows <- rbindlist(result)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(out[[1L]],width=11.69,height=8.27,family="sans",bg="white")
for (w in c(60L,5L)) {
  grid.newpage()
  grid.text(paste0("Supplementary Table 10. Characteristics of Eligible Pairs - ",w,"-Minute Analysis"),
            x=unit(.45,"in"),y=unit(7.72,"in"),just=c("left","centre"),
            gp=gpar(fontsize=11,fontface="bold"))
  xs <- c(.48,4.6,6.1,7.6,9.1,10.35)
  for (j in seq_along(xs)) grid.text(c("Characteristic",unname(pretty[cohorts]))[[j]],
       x=unit(xs[[j]],"in"),y=unit(7.18,"in"),just=c("left","centre"),
       gp=gpar(fontsize=8,fontface="bold"))
  grid.lines(x=unit(c(.45,11.3),"in"),y=unit(c(6.91,6.91),"in"))
  labels <- result[[paste(cohorts[[1L]],w)]]$characteristic
  for (i in seq_along(labels)) {
    y <- 6.54-(i-1)*.43
    grid.text(labels[[i]],x=unit(xs[[1L]],"in"),y=unit(y,"in"),just=c("left","centre"),gp=gpar(fontsize=8))
    for (j in seq_along(cohorts)) {
      val <- result[[paste(cohorts[[j]],w)]]$display_value[[i]]
      grid.text(val,x=unit(xs[[j+1L]],"in"),y=unit(y,"in"),
                just=c("left","centre"),gp=gpar(fontsize=7.7))
    }
  }
  grid.lines(x=unit(c(.45,11.3),"in"),y=unit(c(2.27,2.27),"in"))
}
dev.off()
fwrite(rows,out[[2L]],sep="\t")
fwrite(data.table(table_id="Supplementary Table 10",cohort_n=5L,analysis_n=10L,
                  row_n=nrow(rows)),out[[3L]],sep="\t")
cat("PAIR_CHARACTERISTICS_TABLE_CANDIDATE_PASS rows=",nrow(rows),"\n",sep="")
