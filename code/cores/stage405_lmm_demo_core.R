s405l_stop_if <- function(condition, code) {
  if (isTRUE(condition)) stop(code, call. = FALSE)
  invisible(TRUE)
}

s405l_expected_gate_status <- function() {
  "S405F_UNIFIED_FORMAL_IMPUTATION_INPUT_GATE_PASS_FOR_LMM_NOT_ADOPTED"
}

s405l_required_base_fields <- function() {
  c(
    "bias_spo2_minus_sao2", "paired_mean_saturation_c90",
    "sao2_anchor_icu_hours", "cluster_stay", "cluster_patient"
  )
}

s405l_formula <- function(cohort, predictors) {
  s405m_fixed_formula(cohort, predictors)
}

s405l_fixed_formula <- function(predictors) {
  stats::as.formula(paste0(
    "bias_spo2_minus_sao2 ~ ",
    paste(c("paired_mean_saturation_c90", predictors), collapse = " + ")
  ))
}

s405l_prepare_model_data <- function(d) {
  out <- data.table::copy(d)
  out$cluster_patient <- factor(out$cluster_patient)
  out$cluster_stay <- factor(out$cluster_stay)
  out$sex_class <- stats::relevel(out$sex_class, ref = "male")
  out
}

s405l_fit_attempt <- function(formula, d, optimizer) {
  fit <- tryCatch(
    lme4::lmer(
      formula, data = d, REML = TRUE,
      control = lme4::lmerControl(
        optimizer = optimizer,
        optCtrl = list(maxfun = 100000L),
        calc.derivs = TRUE
      )
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      fit = fit, ok = FALSE, optimizer = optimizer, opt_code = NA_integer_,
      convergence_message = conditionMessage(fit), hessian_available = FALSE,
      hessian_positive_definite = FALSE
    ))
  }
  messages <- fit@optinfo$conv$lme4$messages
  opt_code <- fit@optinfo$conv$opt
  if (is.null(opt_code)) opt_code <- 0L
  hessian <- fit@optinfo$derivs$Hessian
  hessian_available <- is.matrix(hessian) && length(hessian) > 0L && all(is.finite(hessian))
  hessian_positive_definite <- FALSE
  if (hessian_available) {
    eigenvalues <- tryCatch(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values, error = function(e) numeric())
    hessian_positive_definite <- length(eigenvalues) > 0L && all(is.finite(eigenvalues)) && min(eigenvalues) > 0
  }
  ok <- identical(as.integer(opt_code), 0L) && is.null(messages) && hessian_available && hessian_positive_definite
  list(
    fit = fit, ok = ok, optimizer = optimizer, opt_code = as.integer(opt_code),
    convergence_message = paste(messages, collapse = " | "),
    hessian_available = hessian_available,
    hessian_positive_definite = hessian_positive_definite
  )
}

s405l_design_diagnostics <- function(d, predictors) {
  fixed_formula <- s405l_fixed_formula(predictors)
  mm <- Matrix::sparse.model.matrix(fixed_formula, data = d)
  rank <- as.integer(Matrix::rankMatrix(mm, method = "qr")[[1L]])
  data.table::data.table(
    design_column_n = ncol(mm), design_rank = rank,
    full_rank = rank == ncol(mm)
  )
}

s405l_gvif_diagnostics <- function(fit) {
  raw <- suppressMessages(suppressWarnings(car::vif(fit)))
  if (is.matrix(raw)) {
    out <- data.table::data.table(
      term = rownames(raw), gvif = as.numeric(raw[, "GVIF"]),
      df = as.numeric(raw[, "Df"])
    )
    out[, adjusted_gvif := gvif^(1 / (2 * df))]
  } else {
    out <- data.table::data.table(
      term = names(raw), gvif = as.numeric(raw), df = 1,
      adjusted_gvif = sqrt(as.numeric(raw))
    )
  }
  out
}

s405l_residual_diagnostics <- function(fit) {
  raw <- stats::residuals(fit)
  fitted <- stats::fitted(fit)
  sigma <- stats::sigma(fit)
  standardized <- raw / sigma
  probs <- seq(0.005, 0.995, length.out = 199L)
  observed_q <- as.numeric(stats::quantile(standardized, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
  theoretical_q <- stats::qnorm(probs)
  data.table::data.table(
    residual_n = length(raw), all_residuals_finite = all(is.finite(raw)),
    all_fitted_finite = all(is.finite(fitted)), residual_mean = mean(raw),
    residual_sd = stats::sd(raw), standardized_residual_abs_max = max(abs(standardized)),
    residual_fitted_abs_correlation = suppressWarnings(stats::cor(abs(standardized), fitted)),
    qq_quantile_correlation = suppressWarnings(stats::cor(observed_q, theoretical_q))
  )
}

s405l_fit_primary_one <- function(d, cohort, predictors) {
  d <- s405l_prepare_model_data(d)
  formula <- s405l_formula(cohort, predictors)
  first <- s405l_fit_attempt(formula, d, "bobyqa")
  second <- NULL
  chosen <- first
  if (!first$ok) {
    second <- s405l_fit_attempt(formula, d, "nloptwrap")
    chosen <- second
  }
  fit <- chosen$fit
  s405l_stop_if(inherits(fit, "error"), "S405L_PRIMARY_LMM_FIT_FAILURE")
  terms <- names(lme4::fixef(fit))
  variance <- diag(as.matrix(stats::vcov(fit)))[terms]
  vc <- data.table::as.data.table(as.data.frame(lme4::VarCorr(fit)))
  random <- vc[grp != "Residual" & var1 == "(Intercept)" & !is.na(vcov)]
  design <- s405l_design_diagnostics(d, predictors)
  gvif <- s405l_gvif_diagnostics(fit)
  residual <- s405l_residual_diagnostics(fit)
  list(
    formula = formula, terms = terms, coefficients = lme4::fixef(fit)[terms], variances = variance,
    status = data.table::data.table(
      optimizer = chosen$optimizer, optimizer_retry_used = !is.null(second),
      first_optimizer_ok = first$ok, second_optimizer_ok = if (is.null(second)) NA else second$ok,
      opt_code = chosen$opt_code, converged = chosen$ok,
      convergence_message = chosen$convergence_message,
      hessian_available = chosen$hessian_available,
      hessian_positive_definite = chosen$hessian_positive_definite,
      singular = lme4::isSingular(fit, tol = 1e-04), singularity_tolerance = 1e-04,
      random_component_n = nrow(random),
      random_variance_min = if (nrow(random)) min(random$vcov) else NA_real_,
      random_variance_max = if (nrow(random)) max(random$vcov) else NA_real_,
      random_variance_floor_pass = nrow(random) > 0L && all(random$vcov >= 1e-06)
    ),
    random_effects = vc, design = design, gvif = gvif, residual = residual
  )
}

s405l_fit_balanced_one <- function(d, predictors) {
  d <- s405l_prepare_model_data(d)
  formula <- s405l_fixed_formula(predictors)
  patient_n <- table(d$cluster_patient)
  d$patient_balance_weight <- 1 / as.numeric(patient_n[as.character(d$cluster_patient)])
  patient_weight_total <- rowsum(d$patient_balance_weight, d$cluster_patient, reorder = FALSE)
  fit <- stats::lm(formula, data = d, weights = patient_balance_weight)
  covariance <- sandwich::vcovCL(fit, cluster = d$cluster_patient, type = "HC1")
  terms <- names(stats::coef(fit))
  list(
    formula = formula, terms = terms, coefficients = stats::coef(fit)[terms],
    variances = diag(covariance)[terms],
    weight_status = data.table::data.table(
      patient_n = length(patient_n), patient_total_weight_min = min(patient_weight_total),
      patient_total_weight_max = max(patient_weight_total),
      patient_total_weight_max_abs_deviation = max(abs(patient_weight_total - 1)),
      all_patient_total_weights_equal_one = max(abs(patient_weight_total - 1)) <= 1e-10
    )
  )
}

s405l_pool_stream_results <- function(results, patient_n) {
  primary_terms <- results[[1L]]$primary$terms
  balanced_terms <- results[[1L]]$balanced$terms
  s405l_stop_if(!all(vapply(results, function(x) identical(x$primary$terms, primary_terms), logical(1L))), "S405L_PRIMARY_TERM_SET_CHANGED")
  s405l_stop_if(!all(vapply(results, function(x) identical(x$balanced$terms, balanced_terms), logical(1L))), "S405L_BALANCED_TERM_SET_CHANGED")
  s405l_stop_if(!identical(primary_terms, balanced_terms), "S405L_PRIMARY_BALANCED_TERM_SET_DIFFERS")
  pcoef <- do.call(cbind, lapply(results, function(x) x$primary$coefficients[primary_terms]))
  pvar <- do.call(cbind, lapply(results, function(x) x$primary$variances[primary_terms]))
  bcoef <- do.call(cbind, lapply(results, function(x) x$balanced$coefficients[balanced_terms]))
  bvar <- do.call(cbind, lapply(results, function(x) x$balanced$variances[balanced_terms]))
  complete_df <- max(1, patient_n - length(primary_terms))
  primary <- s405m_pool_estimates(pcoef, pvar, complete_df = complete_df)
  balanced <- s405m_pool_estimates(bcoef, bvar, complete_df = complete_df)
  list(
    primary = primary, balanced = balanced,
    comparison = s405m_compare_primary_balanced(primary, balanced),
    primary_coefficients = pcoef, primary_variances = pvar,
    balanced_coefficients = bcoef, balanced_variances = bvar
  )
}

s405l_reporting_scale <- function(term) {
  multiplier <- ifelse(
    term %in% c("pao2_candidate_value", "paco2_candidate_value"), 10,
    ifelse(term == "ph_candidate_value", 0.10, 1)
  )
  label <- ifelse(
    term == "paired_mean_saturation_c90", "Paired mean saturation, per 1 percentage point",
    ifelse(term == "age_value_per_10y", "Age, per 10 years",
    ifelse(term == "age_topcoded_indicator", "Source age top-code indicator",
    ifelse(term == "sex_classfemale", "Female vs male",
    ifelse(term == "cci_value", "Charlson Comorbidity Index, per 1 point",
    ifelse(term == "sofa_total", "SOFA, per 1 point",
    ifelse(term == "pao2_candidate_value", "PaO2, per 10 mm Hg",
    ifelse(term == "paco2_candidate_value", "PaCO2, per 10 mm Hg",
    ifelse(term == "ph_candidate_value", "pH, per 0.10",
    ifelse(term == "lactate_candidate_value", "Lactate, per 1 mmol/L",
    ifelse(term == "hb_candidate_value", "Hemoglobin, per 1 g/dL",
    ifelse(grepl("^age_interval", term), sub("^age_interval", "Age category ", term), term))))))))))))
  data.table::data.table(term = term, reporting_label = label, reporting_multiplier = multiplier)
}

s405l_apply_reporting_scale <- function(pooled) {
  scale <- s405l_reporting_scale(pooled$term)
  out <- merge(data.table::copy(pooled), scale, by = "term", sort = FALSE)
  for (field in intersect(c("estimate", "std_error", "ci_lower", "ci_upper", "monte_carlo_error_estimate"), names(out))) {
    out[, (field) := get(field) * reporting_multiplier]
  }
  out[]
}
