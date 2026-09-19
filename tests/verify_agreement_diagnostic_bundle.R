#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: DIAGNOSTIC_MANIFEST.tsv PROPORTIONAL_ROWS.tsv HETEROSCEDASTIC_ROWS.tsv",call.=FALSE)
m<-fread(args[[1L]])
p<-fread(args[[2L]],colClasses="character")
h<-fread(args[[3L]],colClasses="character")
pretty<-c(MIMIC="MIMIC",Amsterdam="AmsterdamUMCdb",eICU="eICU",
          SICDB="SICdb",Lianyungang="Lianyungang")
if(nrow(m)!=10L || nrow(p)!=10L || nrow(h)!=10L) stop("DIAGNOSTIC_BUNDLE_QA_SHAPE_FAIL",call.=FALSE)
checks<-0L
for(i in seq_len(nrow(m))) {
  co<-m$cohort[[i]]
  w<-m$window_minutes[[i]]
  pp<-fread(file.path(m$diagnostic_dir[[i]],"agreement_proportional_bias.tsv"))
  hh<-fread(file.path(m$diagnostic_dir[[i]],"agreement_heteroscedasticity.tsv"))
  rp<-p[cohort==pretty[[co]] & window==paste0(w," min")]
  rh<-h[cohort==pretty[[co]] & window==paste0(w," min")]
  if(nrow(rp)!=1L || nrow(rh)!=1L ||
     rp$slope_with_95ci!=sprintf("%.2f (%.2f to %.2f)",pp$slope,pp$ci_lower,pp$ci_upper) ||
     rp$model_estimated_bias_at_90!=sprintf("%.2f",pp$model_estimated_bias_at_90) ||
     rh$variance_parameter_with_95ci!=sprintf("%.3f (%.3f to %.3f)",
        hh$variance_function_parameter,hh$ci_lower,hh$ci_upper) ||
     rh$residual_sd_at_90!=sprintf("%.2f",hh$residual_sd_at_90) ||
     rh$likelihood_ratio!=sprintf("%.1f",hh$likelihood_ratio) ||
     rh$p_value!=if(hh$p_value<.001) "<.001" else sprintf("%.3f",hh$p_value))
    stop(paste("DIAGNOSTIC_BUNDLE_QA_VALUE_FAIL",co,w),call.=FALSE)
  checks<-checks+6L
}
cat("DIAGNOSTIC_BUNDLE_INDEPENDENT_QA_PASS checks=",checks,"\n",sep="")
