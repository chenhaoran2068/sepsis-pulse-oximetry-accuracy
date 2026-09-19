#!/usr/bin/env Rscript
# Independent invented-data check of model-estimated Bias at paired mean 90%.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: AGREEMENT_CORE.R", call. = FALSE)
source(args[[1L]])
suppressPackageStartupMessages(library(data.table))
set.seed(91821)
z <- CJ(patient = sprintf("p%02d", 1:30), centred_mean = seq(-8, 8, length.out = 15))
z[, patient := factor(patient)]
patient_effect <- rnorm(30, sd = 0.5)
z[, bias := 0.75 + 0.12 * centred_mean + patient_effect[as.integer(patient)] +
     rnorm(.N, sd = 0.4)]
fit <- fit_proportional_bias(z, nested = FALSE)
reference_intercept <- unname(lme4::fixef(fit$model)[["(Intercept)"]])
reference_at_90 <- unname(predict(fit$model,
                                  newdata = data.frame(centred_mean = 0, patient = z$patient[[1L]]),
                                  re.form = NA))
if (abs(fit$intercept_at_90 - reference_intercept) > 1e-12 ||
    abs(fit$intercept_at_90 - reference_at_90) > 1e-12 ||
    abs(fit$intercept_at_90 - 0.75) > 0.25 || abs(fit$slope - 0.12) > 0.02) {
  stop("PROPORTIONAL_INTERCEPT_OR_SLOPE_WRONG", call. = FALSE)
}
cat("PROPORTIONAL_INTERCEPT_INDEPENDENT_QA_PASS intercept=", fit$intercept_at_90,
    " slope=", fit$slope, "\n", sep = "")
