#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: LMM60_MANIFEST.tsv LMM5_MANIFEST.tsv ROWS.tsv",call.=FALSE)
d <- fread(args[[3L]])
pretty <- c(MIMIC="MIMIC",Amsterdam="AmsterdamUMCdb",eICU="eICU",
            SICDB="SICdb",Lianyungang="Lianyungang")
if(nrow(d)<100L || anyDuplicated(d[,.(cohort,window_minutes,feature,term)]))
  stop("LMM_CHARACTERISTICS_QA_SHAPE_FAIL",call.=FALSE)
checks <- 0L
for(k in 1:2) {
  w <- c(60L,5L)[[k]]
  m <- fread(args[[k]])
  for(i in seq_len(nrow(m))) {
    co <- m$cohort[[i]]
    panel <- d[cohort==pretty[[co]] & window_minutes==w]
    agg <- dirname(m$lmm_tsv[[i]])
    p <- fread(m$lmm_tsv[[i]])
    cvr <- fread(file.path(agg,"variable_coverage.tsv"))
    run <- fread(file.path(agg,"run_summary.tsv"))
    expected_terms <- setdiff(p$term,"(Intercept)")
    populated_terms <- panel$term[!is.na(panel$term) & nzchar(panel$term)]
    if(!setequal(populated_terms,expected_terms) ||
       anyDuplicated(populated_terms))
      stop(paste("LMM_CHARACTERISTICS_QA_TERM_SET_FAIL",co,w),call.=FALSE)
    for(j in seq_len(nrow(panel))) {
      row <- panel[j]
      if(!is.na(row$term) && nzchar(row$term)) {
        q <- p[match(row$term,p$term)]
        mul <- q$reporting_multiplier
        if(abs(row$estimate-q$estimate*mul)>1e-10 ||
           abs(row$ci_lower-q$ci_lower*mul)>1e-10 ||
           abs(row$ci_upper-q$ci_upper*mul)>1e-10 ||
           abs(row$p_value-q$p_value)>1e-10)
          stop("LMM_CHARACTERISTICS_QA_ESTIMATE_FAIL",call.=FALSE)
        checks <- checks+4L
      }
      if(!is.na(row$observed_n)) {
        key <- row$feature
        if(key=="age") key <- if("age_interval" %in% cvr$feature)
          "age_interval" else "age_value_per_10y"
        if(key=="paired_mean_saturation_c90") {
          n <- run$pair_n
          denominator <- n
        } else {
          q <- cvr[feature==key]
          if(nrow(q)!=1L) stop("LMM_CHARACTERISTICS_QA_COVERAGE_KEY_FAIL",call.=FALSE)
          n <- q$observed_n
          denominator <- q$denominator_n
        }
        if(row$observed_n!=n || row$denominator_n!=denominator ||
           abs(row$observed_pct-100*n/denominator)>1e-10)
          stop("LMM_CHARACTERISTICS_QA_COVERAGE_FAIL",call.=FALSE)
        checks <- checks+3L
      }
    }
  }
}
cat("LMM_CHARACTERISTICS_INDEPENDENT_QA_PASS checks=",checks,"\n",sep="")
