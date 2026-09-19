#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: LMM60_MANIFEST.tsv LMM5_MANIFEST.tsv ROWS.tsv",call.=FALSE)
d <- fread(args[[3L]])
if(nrow(d)!=10L || anyDuplicated(d[,.(cohort,window_minutes)]))
  stop("LMM_DIAGNOSTICS_QA_SHAPE_FAIL",call.=FALSE)
pretty <- c(MIMIC="MIMIC",Amsterdam="AmsterdamUMCdb",eICU="eICU",
            SICDB="SICdb",Lianyungang="Lianyungang")
checks <- 0L
for(k in 1:2) {
  w <- c(60L,5L)[[k]]
  m <- fread(args[[k]])
  for(i in seq_len(nrow(m))) {
    cohort_value <- m$cohort[[i]]
    row <- d[cohort==pretty[[cohort_value]] & window_minutes==w]
    agg <- dirname(m$lmm_tsv[[i]])
    fits <- fread(file.path(agg,"primary_fit_status.tsv"))
    gvif <- fread(file.path(agg,"primary_gvif.tsv"))
    residual <- fread(file.path(agg,"primary_residual_status.tsv"))
    run <- fread(file.path(agg,"run_summary.tsv"))
    if(nrow(row)!=1L || row$total_imputations!=run$imputations ||
       row$successfully_fitted_lmms!=sum(fits$converged) ||
       row$singular_fits!=sum(fits$singular) ||
       abs(row$largest_adjusted_gvif-max(gvif$adjusted_gvif))>1e-10 ||
       abs(row$residual_fitted_abs_correlation_min-
           min(residual$residual_fitted_abs_correlation))>1e-10 ||
       abs(row$residual_fitted_abs_correlation_max-
           max(residual$residual_fitted_abs_correlation))>1e-10)
      stop(paste("LMM_DIAGNOSTICS_QA_VALUE_FAIL",cohort_value,w),call.=FALSE)
    group_variances <- list()
    for(j in seq_len(run$imputations)) {
      f <- readRDS(file.path(dirname(agg),"checkpoints",
                             sprintf("fit_%02d.rds",j)))
      re <- as.data.table(f$primary$random_effects)
      group_variances[[j]] <- re[var1=="(Intercept)" & is.na(var2),.(grp,vcov)]
    }
    pooled <- rbindlist(group_variances)[,.(lo=min(vcov),hi=max(vcov)),by=grp]
    for(j in seq_len(nrow(pooled))) {
      token <- sprintf("%s: %.3f to %.3f",pooled$grp[[j]],
                       pooled$lo[[j]],pooled$hi[[j]])
      if(!grepl(token,row$random_effect_variance_range,fixed=TRUE))
        stop(paste("LMM_DIAGNOSTICS_QA_VARIANCE_FAIL",cohort_value,w),call.=FALSE)
      checks <- checks+1L
    }
    checks <- checks+6L
  }
}
cat("LMM_DIAGNOSTICS_INDEPENDENT_QA_PASS checks=",checks,"\n",sep="")
