#!/usr/bin/env Rscript
# Provenance gate: all selected outputs must stem from the same analysis-input files.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(digest))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: FIVE_COHORT_RUN_MAP.tsv",call.=FALSE)
m <- fread(args[[1L]])
cohorts <- c("MIMIC","Amsterdam","eICU","SICDB","Lianyungang")
fields <- c("cohort","analysis_dir","agreement_60_dir","agreement_5_dir",
            "threshold_60_dir","threshold_5_dir","lmm_60_dir","lmm_5_dir","rf_60_dir")
if(!identical(names(m),fields) || nrow(m)!=5L ||
   !setequal(m$cohort,cohorts) || anyDuplicated(m$cohort))
  stop("DISPLAY_RUN_MAP_INVALID",call.=FALSE)
resolve <- function(path) {
  if(grepl("^[A-Za-z]:[/\\]|^/",path)) return(path)
  file.path(dirname(normalizePath(args[[1L]],winslash="/",mustWork=TRUE)),path)
}
hash <- function(path) digest(path,file=TRUE,algo="sha256")
checks <- 0L
for(i in seq_len(nrow(m))) {
  co <- m$cohort[[i]]
  analysis <- resolve(m$analysis_dir[[i]])
  for(w in c(60L,5L)) {
    analytic_path <- file.path(analysis,paste0("analytic_pairs_",w,".rds"))
    matched_path <- file.path(analysis,paste0("matched_pairs_",w,".rds"))
    lmm_path <- file.path(analysis,paste0("lmm_input_",w,".rds"))
    agreement_dir <- resolve(m[[paste0("agreement_",w,"_dir")]][[i]])
    threshold_dir <- resolve(m[[paste0("threshold_",w,"_dir")]][[i]])
    lmm_dir <- resolve(m[[paste0("lmm_",w,"_dir")]][[i]])
    agreement <- fread(file.path(agreement_dir,"run_receipt.tsv"))
    threshold <- fread(file.path(threshold_dir,"run_receipt.tsv"))
    lmm <- fread(file.path(lmm_dir,"aggregate","run_summary.tsv"))
    if(!file.exists(analytic_path) || !file.exists(matched_path) ||
       !file.exists(lmm_path) ||
       nrow(agreement)!=1L || nrow(threshold)!=1L || nrow(lmm)!=1L ||
       agreement$cohort[[1L]]!=co || threshold$cohort[[1L]]!=co ||
       lmm$cohort[[1L]]!=co ||
       agreement$window_minutes[[1L]]!=w ||
       threshold$window_minutes[[1L]]!=w || lmm$window_minutes[[1L]]!=w ||
       agreement$input_sha256[[1L]]!=hash(matched_path) ||
       threshold$full_sha256[[1L]]!=hash(matched_path) ||
       lmm$input_sha256[[1L]]!=hash(lmm_path))
      stop(paste("DISPLAY_RUN_MAP_PROVENANCE_MISMATCH",co,w),call.=FALSE)
    checks <- checks+9L
  }
  rf <- readRDS(file.path(resolve(m$rf_60_dir[[i]]),"run_identity.rds"))
  candidate <- file.path(analysis,"rf_candidate_60.rds")
  if(!file.exists(candidate) || rf$cohort!=co ||
     rf$candidate_sha256!=hash(candidate))
    stop(paste("DISPLAY_RUN_MAP_RF_PROVENANCE_MISMATCH",co),call.=FALSE)
  checks <- checks+2L
}
cat("DISPLAY_RUN_MAP_PROVENANCE_PASS checks=",checks,"\n",sep="")
