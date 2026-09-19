#!/usr/bin/env Rscript
# Technical renderer for Supplementary Table 20. Does not alter formal model fits.
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
features <- c("paired_mean_saturation_c90","age","age_topcoded_indicator","sex_class",
              "cci_value","sofa_total","pao2_candidate_value","paco2_candidate_value",
              "ph_candidate_value","lactate_candidate_value","hb_candidate_value")
labels <- c("Paired mean saturation","Age","Source top-coded age indicator",
            "Sex","Charlson Comorbidity Index","SOFA","PaO2","PaCO2","pH",
            "Lactate","Hemoglobin")
map_term <- function(feature,pooled) {
  switch(feature,
    paired_mean_saturation_c90="paired_mean_saturation_c90",
    age=pooled$term[grepl("^age_value_per_10y$|^age_interval",pooled$term)],
    age_topcoded_indicator="age_topcoded_indicator",
    sex_class="sex_classfemale",feature)
}
rows <- list()
for(k in seq_along(c(60L,5L))) {
  w <- c(60L,5L)[[k]]
  m <- fread(args[[k]])
  if(!identical(names(m),c("cohort","lmm_tsv")) || nrow(m)!=5L ||
     !setequal(m$cohort,cohorts) || anyDuplicated(m$cohort))
    stop("LMM_CHARACTERISTICS_MANIFEST_INVALID",call.=FALSE)
  for(cohort_value in cohorts) {
    pooled_path <- m[cohort==cohort_value]$lmm_tsv
    agg <- dirname(pooled_path)
    coverage_path <- file.path(agg,"variable_coverage.tsv")
    run_path <- file.path(agg,"run_summary.tsv")
    if(length(pooled_path)!=1L || !all(file.exists(c(pooled_path,coverage_path,run_path))))
      stop("LMM_CHARACTERISTICS_SOURCE_MISSING",call.=FALSE)
    p <- fread(pooled_path)
    cvr <- fread(coverage_path)
    run <- fread(run_path)
    if(anyDuplicated(p$term) || anyDuplicated(cvr$feature) || nrow(run)!=1L ||
       run$cohort[[1L]]!=cohort_value || run$window_minutes[[1L]]!=w ||
       !all(c("term","estimate","ci_lower","ci_upper","p_value",
              "reporting_label","reporting_multiplier") %in% names(p)) ||
       !all(c("feature","observed_n","denominator_n","core_model_eligible") %in% names(cvr)))
      stop("LMM_CHARACTERISTICS_SOURCE_INVALID",call.=FALSE)
    panel <- list()
    for(i in seq_along(features)) {
      feature <- features[[i]]
      coverage_key <- if(feature=="age") {
        if("age_interval" %in% cvr$feature) "age_interval" else "age_value_per_10y"
      } else feature
      if(feature=="paired_mean_saturation_c90") {
        n <- run$pair_n[[1L]]
        denominator <- n
        eligible <- TRUE
      } else if(coverage_key %in% cvr$feature) {
        rr <- cvr[feature==coverage_key]
        n <- rr$observed_n[[1L]]
        denominator <- rr$denominator_n[[1L]]
        eligible <- isTRUE(rr$core_model_eligible[[1L]])
      } else if(feature=="age_topcoded_indicator") {
        next
      } else stop(paste("LMM_CHARACTERISTICS_COVERAGE_MISSING",feature),call.=FALSE)
      terms <- intersect(map_term(feature,p),p$term)
      if(eligible && length(terms)==0L && feature!="age")
        stop(paste("LMM_CHARACTERISTICS_MODEL_TERM_MISSING",feature),call.=FALSE)
      if(!eligible && length(terms)>0L)
        stop(paste("LMM_CHARACTERISTICS_INELIGIBLE_TERM_PRESENT",feature),call.=FALSE)
      if(length(terms)==0L) terms <- NA_character_
      for(j in seq_along(terms)) {
        term <- terms[[j]]
        fitted <- if(is.na(term)) NULL else p[match(terms[[j]],p$term)]
        if(!is.null(fitted) && nrow(fitted)!=1L)
          stop("LMM_CHARACTERISTICS_TERM_AMBIGUOUS",call.=FALSE)
        est <- if(is.null(fitted)) NA_real_ else fitted$estimate*fitted$reporting_multiplier
        lo <- if(is.null(fitted)) NA_real_ else fitted$ci_lower*fitted$reporting_multiplier
        hi <- if(is.null(fitted)) NA_real_ else fitted$ci_upper*fitted$reporting_multiplier
        pval <- if(is.null(fitted)) NA_real_ else fitted$p_value
        if(!is.null(fitted) && (any(!is.finite(c(est,lo,hi,pval))) || lo>est || est>hi))
          stop("LMM_CHARACTERISTICS_FIT_VALUE_INVALID",call.=FALSE)
        panel[[paste(feature,j)]] <- data.table(
          cohort=unname(pretty[[cohort_value]]),window_minutes=w,
          feature=feature,term=term,
          label=if(!is.null(fitted)) fitted$reporting_label else labels[[i]],
          observed_n=if(j==1L) n else NA_integer_,
          denominator_n=if(j==1L) denominator else NA_integer_,
          observed_pct=if(j==1L) 100*n/denominator else NA_real_,
          model_included=eligible,estimate=est,ci_lower=lo,ci_upper=hi,p_value=pval)
      }
    }
    rows[[paste(cohort_value,w)]] <- rbindlist(panel)
  }
}
d <- rbindlist(rows,fill=TRUE)
pdf_device <- if(capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(out[[1L]],width=11.69,height=8.27,family="sans",bg="white")
for(w in c(60L,5L)) for(cohort_value in cohorts) {
  panel <- rows[[paste(cohort_value,w)]]
  if(nrow(panel)>27L) stop("LMM_CHARACTERISTICS_TOO_MANY_ROWS",call.=FALSE)
  grid.newpage()
  grid.text(paste0("Supplementary Table 20. Observed Coverage and LMM Associations - ",
                   unname(pretty[[cohort_value]])," / ",w," min"),
            x=unit(.43,"in"),y=unit(7.72,"in"),just=c("left","centre"),
            gp=gpar(fontsize=10.2,fontface="bold"))
  xs <- c(.45,5.55,7.7,10.3)
  headers <- c("Variable, comparison, or unit","Observed, n/N (%)",
               "Coefficient (95% CI)","P value")
  for(j in seq_along(xs)) grid.text(headers[[j]],x=unit(xs[[j]],"in"),
       y=unit(7.2,"in"),just=c("left","centre"),gp=gpar(fontsize=7.8,fontface="bold"))
  grid.lines(x=unit(c(.43,11.3),"in"),y=unit(c(6.96,6.96),"in"))
  row_step <- min(.28,5.35/max(1L,nrow(panel)-1L))
  ys <- seq(6.68,by=-row_step,length.out=nrow(panel))
  for(i in seq_len(nrow(panel))) {
    r <- panel[i]
    coverage <- if(is.na(r$observed_n)) "-" else
      sprintf("%d/%d (%.1f)",r$observed_n,r$denominator_n,r$observed_pct)
    effect <- if(is.na(r$estimate)) "-" else
      sprintf("%+.2f (%+.2f to %+.2f)",r$estimate,r$ci_lower,r$ci_upper)
    pval <- if(is.na(r$p_value)) "-" else
      if(r$p_value<.001) "<.001" else sprintf("%.3f",r$p_value)
    vals <- c(r$label,coverage,effect,pval)
    for(j in seq_along(xs)) grid.text(vals[[j]],x=unit(xs[[j]],"in"),
         y=unit(ys[[i]],"in"),just=c("left","centre"),gp=gpar(fontsize=7.1))
  }
  grid.lines(x=unit(c(.43,11.3),"in"),
             y=unit(c(min(ys)-.19,min(ys)-.19),"in"))
}
dev.off()
fwrite(d,out[[2L]],sep="\t")
fwrite(data.table(table_id="Supplementary Table 20",cohort_n=5L,analysis_n=10L,
                  row_n=nrow(d)),out[[3L]],sep="\t")
cat("LMM_CHARACTERISTICS_TABLE_CANDIDATE_PASS rows=",nrow(d),"\n",sep="")
