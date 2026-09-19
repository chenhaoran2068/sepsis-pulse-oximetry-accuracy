#!/usr/bin/env Rscript
# Technical renderer for Supplementary Table 21 from formal model aggregates/checkpoints.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop("Usage: LMM60_MANIFEST.tsv LMM5_MANIFEST.tsv NEW.pdf NEW_ROWS.tsv NEW_RECEIPT.tsv",call.=FALSE)
cohorts <- c("MIMIC","Amsterdam","eICU","SICDB","Lianyungang")
pretty <- c(MIMIC="MIMIC",Amsterdam="AmsterdamUMCdb",eICU="eICU",
            SICDB="SICdb",Lianyungang="Lianyungang")
out <- args[3:5]
if(any(file.exists(out))) stop("OUTPUT_EXISTS_REFUSE_OVERWRITE",call.=FALSE)
if(any(!dir.exists(dirname(out)))) stop("OUTPUT_PARENT_MISSING",call.=FALSE)
rows <- list()
for(k in seq_along(c(60L,5L))) {
  w <- c(60L,5L)[[k]]
  m <- fread(args[[k]])
  if(!identical(names(m),c("cohort","lmm_tsv")) || nrow(m)!=5L ||
     !setequal(m$cohort,cohorts) || anyDuplicated(m$cohort))
    stop("LMM_DIAGNOSTIC_MANIFEST_INVALID",call.=FALSE)
  for(cohort_value in cohorts) {
    pooled <- m[cohort==cohort_value]$lmm_tsv
    agg <- dirname(pooled)
    run_dir <- dirname(agg)
    paths <- file.path(agg,c("primary_fit_status.tsv","primary_gvif.tsv",
                             "primary_residual_status.tsv","run_summary.tsv"))
    if(length(pooled)!=1L || !all(file.exists(paths)))
      stop("LMM_DIAGNOSTIC_SOURCE_MISSING",call.=FALSE)
    fits <- fread(paths[[1L]])
    gvif <- fread(paths[[2L]])
    residual <- fread(paths[[3L]])
    run <- fread(paths[[4L]])
    expected <- as.integer(run$imputations[[1L]])
    if(nrow(run)!=1L || run$cohort[[1L]]!=cohort_value ||
       run$window_minutes[[1L]]!=w || expected<1L ||
       nrow(fits)!=expected || nrow(residual)!=expected ||
       !setequal(fits$imputation,seq_len(expected)) ||
       !setequal(residual$imputation,seq_len(expected)) ||
       !setequal(gvif$imputation,seq_len(expected)) ||
       anyNA(fits[,.(converged,singular)]) ||
       any(!is.finite(gvif$adjusted_gvif)) ||
       any(!is.finite(residual$residual_fitted_abs_correlation)))
      stop("LMM_DIAGNOSTIC_SOURCE_INVALID",call.=FALSE)
    random_effects <- vector("list",expected)
    for(i in seq_len(expected)) {
      path <- file.path(run_dir,"checkpoints",sprintf("fit_%02d.rds",i))
      if(!file.exists(path)) stop("LMM_DIAGNOSTIC_CHECKPOINT_MISSING",call.=FALSE)
      x <- readRDS(path)
      re <- as.data.table(x$primary$random_effects)
      re <- re[var1=="(Intercept)" & is.na(var2),.(grp,vcov)]
      if(nrow(re)<1L || anyNA(re) || any(re$vcov<0) || anyDuplicated(re$grp))
        stop("LMM_DIAGNOSTIC_RANDOM_EFFECT_INVALID",call.=FALSE)
      re[,imputation:=i]
      random_effects[[i]] <- re
    }
    rv <- rbindlist(random_effects)
    if(any(rv[,uniqueN(grp),by=imputation]$V1 != uniqueN(rv$grp)))
      stop("LMM_DIAGNOSTIC_RANDOM_EFFECT_COMPONENT_MISMATCH",call.=FALSE)
    ranges <- rv[,.(minimum=min(vcov),maximum=max(vcov)),by=grp]
    if(cohort_value=="eICU" && nrow(ranges)!=2L)
      stop("LMM_DIAGNOSTIC_EICU_COMPONENT_MISSING",call.=FALSE)
    range_text <- paste(vapply(seq_len(nrow(ranges)),function(i)
      sprintf("%s: %.3f to %.3f",ranges$grp[[i]],ranges$minimum[[i]],
              ranges$maximum[[i]]),character(1)),collapse="; ")
    cor_min <- min(residual$residual_fitted_abs_correlation)
    cor_max <- max(residual$residual_fitted_abs_correlation)
    rows[[paste(cohort_value,w)]] <- data.table(
      cohort=unname(pretty[[cohort_value]]),window_minutes=w,
      successfully_fitted_lmms=sum(fits$converged),total_imputations=expected,
      singular_fits=sum(fits$singular),largest_adjusted_gvif=max(gvif$adjusted_gvif),
      random_effect_variance_range=range_text,
      residual_fitted_abs_correlation_min=cor_min,
      residual_fitted_abs_correlation_max=cor_max)
  }
}
d <- rbindlist(rows)
pdf_device <- if(capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(out[[1L]],width=11.69,height=8.27,family="sans",bg="white")
grid.newpage()
grid.text("Supplementary Table 21. Multiple-Imputation and LMM Fit Diagnostics",
          x=unit(.4,"in"),y=unit(7.7,"in"),just=c("left","centre"),
          gp=gpar(fontsize=10.5,fontface="bold"))
xs <- c(.45,2.3,3.45,4.55,5.65,7.4,10.25)
headers <- c("Cohort","Analysis","Fitted LMMs","Singular fits","Largest adj. GVIF",
             "Random-effect variance range","Residual correlation")
for(j in seq_along(xs)) grid.text(headers[[j]],x=unit(xs[[j]],"in"),
     y=unit(7.12,"in"),just=c("left","centre"),gp=gpar(fontsize=7.2,fontface="bold"))
grid.lines(x=unit(c(.4,11.3),"in"),y=unit(c(6.84,6.84),"in"))
for(i in seq_len(nrow(d))) {
  y <- 6.47-(i-1)*.54
  vals <- c(d$cohort[[i]],paste0(d$window_minutes[[i]]," min"),
            paste0(d$successfully_fitted_lmms[[i]],"/",d$total_imputations[[i]]),
            paste0(d$singular_fits[[i]],"/",d$total_imputations[[i]]),
            sprintf("%.2f",d$largest_adjusted_gvif[[i]]),
            gsub("cluster_patient","Patient",
                 gsub("cluster_stay:cluster_patient","ICU stay",
                      d$random_effect_variance_range[[i]],fixed=TRUE),fixed=TRUE),
            if(abs(d$residual_fitted_abs_correlation_min[[i]]-
                   d$residual_fitted_abs_correlation_max[[i]])<.005)
              sprintf("%.2f",d$residual_fitted_abs_correlation_min[[i]]) else
              sprintf("%.2f to %.2f",d$residual_fitted_abs_correlation_min[[i]],
                      d$residual_fitted_abs_correlation_max[[i]]))
  for(j in seq_along(xs)) grid.text(vals[[j]],x=unit(xs[[j]],"in"),
       y=unit(y,"in"),just=c("left","centre"),gp=gpar(fontsize=6.7))
}
grid.lines(x=unit(c(.4,11.3),"in"),y=unit(c(.98,.98),"in"))
dev.off()
fwrite(d,out[[2L]],sep="\t")
fwrite(data.table(table_id="Supplementary Table 21",cohort_n=5L,analysis_n=10L,
                  row_n=nrow(d)),out[[3L]],sep="\t")
cat("LMM_DIAGNOSTICS_TABLE_CANDIDATE_PASS rows=",nrow(d),"\n",sep="")
