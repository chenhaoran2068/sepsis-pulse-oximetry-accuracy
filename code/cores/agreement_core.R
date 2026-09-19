suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

required_pair_columns <- c(
  "pair_key_internal", "cohort", "patient_key_internal", "stay_key_internal",
  "spo2_saturation_percent", "sao2_saturation_percent",
  "absolute_lag_minutes", "pair_direction", "pair_window_minutes",
  "contains_possible_transient_artifact"
)

validate_pair_input <- function(x, cohort_name, window_name) {
  miss <- setdiff(required_pair_columns, names(x))
  if (length(miss)) stop("PAIR_INPUT_COLUMNS_MISSING", call. = FALSE)
  if (!nrow(x)) stop("PAIR_INPUT_EMPTY", call. = FALSE)
  if (anyNA(x[, ..required_pair_columns])) stop("PAIR_INPUT_REQUIRED_VALUE_MISSING", call. = FALSE)
  if (uniqueN(x$pair_key_internal) != nrow(x)) stop("PAIR_KEY_NOT_UNIQUE", call. = FALSE)
  if (uniqueN(x$cohort) != 1L) stop("PAIR_COHORT_NOT_UNIQUE", call. = FALSE)
  if (as.character(unique(x$cohort)) != cohort_name) stop("PAIR_COHORT_MISMATCH", call. = FALSE)
  expected_minutes <- if (window_name == "M60") 60 else 5
  if (!all(as.numeric(x$pair_window_minutes) == expected_minutes)) stop("PAIR_WINDOW_MISMATCH", call. = FALSE)
  if (any(as.numeric(x$absolute_lag_minutes) < 0 | as.numeric(x$absolute_lag_minutes) > expected_minutes + 1e-8)) stop("PAIR_LAG_OUTSIDE_WINDOW", call. = FALSE)
  if (any(as.numeric(x$spo2_saturation_percent) < 70 | as.numeric(x$spo2_saturation_percent) > 100)) stop("SPO2_OUTSIDE_ADOPTED_RANGE", call. = FALSE)
  if (any(as.numeric(x$sao2_saturation_percent) < 70 | as.numeric(x$sao2_saturation_percent) > 100)) stop("SAO2_OUTSIDE_ADOPTED_RANGE", call. = FALSE)
  if (any(as.logical(x$contains_possible_transient_artifact), na.rm = TRUE)) stop("FINAL_PAIR_CONTAINS_ARTIFACT", call. = FALSE)
  invisible(TRUE)
}

prepare_agreement_data <- function(x) {
  z <- x[, .(
    pair_key_internal = as.character(pair_key_internal),
    patient = as.character(patient_key_internal),
    stay = as.character(stay_key_internal),
    spo2 = as.numeric(spo2_saturation_percent),
    sao2 = as.numeric(sao2_saturation_percent)
  )]
  z[, bias := spo2 - sao2]
  z[, paired_mean := (spo2 + sao2) / 2]
  z[, centred_mean := paired_mean - 90]
  z[, spo2_stratum := fifelse(
    spo2 < 88, "70%-<88%",
    fifelse(spo2 < 92, "88%-<92%", fifelse(spo2 <= 96, "92%-<=96%", ">96%-100%"))
  )]
  z[, patient_index := match(patient, unique(patient))]
  z[]
}

cluster_summaries <- function(z) {
  setorder(z, patient_index)
  z[, .(
    n = .N,
    mean_y = mean(bias),
    mean_y2 = mean(bias^2),
    within_ss = sum((bias - mean(bias))^2)
  ), by = patient_index]
}

nested_stay_summaries <- function(z) {
  setorder(z, patient_index, stay)
  z[, .(
    n = .N,
    mean_y = mean(bias),
    within_ss = sum((bias - mean(bias))^2)
  ), by = .(patient_index, stay)]
}

nested_reml_components <- function(log_par, ss) {
  tau2 <- exp(log_par[[1L]])
  gamma2 <- exp(log_par[[2L]])
  sigma2 <- exp(log_par[[3L]])
  denom <- sigma2 + ss$n * gamma2
  a <- ss$n / denom
  aggregate_matrix <- cbind(
    a = a,
    b = a * ss$mean_y,
    c = a * ss$mean_y^2,
    logdet_a = (ss$n - 1) * log(sigma2) + log(denom),
    within = ss$within_ss / sigma2,
    n = ss$n
  )
  aggregated <- rowsum(aggregate_matrix, group = ss$patient_index, reorder = FALSE)
  per_patient <- as.data.table(aggregated)
  per_patient[, patient_index := seq_len(.N)]
  one_plus <- 1 + tau2 * per_patient$a
  per_patient[, `:=`(
    info = a / one_plus,
    xvy = b / one_plus,
    logdet = logdet_a + log(one_plus)
  )]
  list(tau2 = tau2, gamma2 = gamma2, sigma2 = sigma2, per_patient = per_patient)
}

nested_reml_objective <- function(log_par, ss, freq) {
  comp <- nested_reml_components(log_par, ss)
  pp <- comp$per_patient
  if (length(freq) != nrow(pp)) return(.Machine$double.xmax / 100)
  info <- sum(freq * pp$info)
  n_total <- sum(freq * pp$n)
  if (!is.finite(info) || info <= 0 || n_total <= 1) return(.Machine$double.xmax / 100)
  mu <- sum(freq * pp$xvy) / info
  centred_b <- pp$b - mu * pp$a
  quad_patient <- pp$within + pp$c - 2 * mu * pp$b + mu^2 * pp$a -
    comp$tau2 * centred_b^2 / (1 + comp$tau2 * pp$a)
  logdet <- sum(freq * pp$logdet)
  quad <- sum(freq * quad_patient)
  0.5 * (logdet + log(info) + quad + (n_total - 1) * log(2 * pi))
}

fit_nested_reml_sufficient <- function(ss, patient_count, freq = rep.int(1, patient_count), start = NULL, require_hessian = FALSE) {
  if (length(freq) != patient_count) stop("NESTED_FREQUENCY_LENGTH_MISMATCH", call. = FALSE)
  if (sum(freq > 0) < 2L || sum(freq[ss$patient_index] * ss$n) < 3L) stop("INSUFFICIENT_NESTED_CLUSTERS", call. = FALSE)
  if (is.null(start)) {
    weighted_n <- freq[ss$patient_index] * ss$n
    total_mean <- sum(weighted_n * ss$mean_y) / sum(weighted_n)
    total_ss <- sum(freq[ss$patient_index] * ss$within_ss) + sum(weighted_n * (ss$mean_y - total_mean)^2)
    total_var <- max(total_ss / max(sum(weighted_n) - 1, 1), 1e-4)
    sigma2_0 <- max(
      sum(freq[ss$patient_index] * ss$within_ss) /
        max(sum(freq[ss$patient_index] * (ss$n - 1)), 1),
      total_var * 0.25, 1e-5
    )
    remaining <- max(total_var - sigma2_0, total_var * 0.2)
    start <- log(c(max(remaining * 0.6, 1e-5), max(remaining * 0.4, 1e-5), sigma2_0))
  }
  opt <- optim(
    par = start, fn = nested_reml_objective, ss = ss, freq = as.numeric(freq),
    method = "L-BFGS-B", lower = log(rep(1e-10, 3L)), upper = log(rep(1e4, 3L)),
    control = list(maxit = 600, factr = 1e3, pgtol = 1e-8)
  )
  if (opt$convergence != 0L || !is.finite(opt$value)) stop("NESTED_REML_OPTIMIZATION_FAILED", call. = FALSE)
  comp <- nested_reml_components(opt$par, ss)
  pp <- comp$per_patient
  info <- sum(freq * pp$info)
  mu <- sum(freq * pp$xvy) / info
  total_sd <- sqrt(comp$tau2 + comp$gamma2 + comp$sigma2)
  hessian_min_eigen <- NA_real_
  if (require_hessian) {
    h <- optimHess(opt$par, nested_reml_objective, ss = ss, freq = as.numeric(freq))
    ev <- eigen((h + t(h)) / 2, symmetric = TRUE, only.values = TRUE)$values
    hessian_min_eigen <- min(ev)
    if (any(!is.finite(ev)) || hessian_min_eigen <= 0) stop("NESTED_REML_HESSIAN_NOT_POSITIVE", call. = FALSE)
  }
  list(
    mean_bias = mu,
    between_patient_variance = comp$tau2,
    between_stay_within_patient_variance = comp$gamma2,
    within_stay_residual_variance = comp$sigma2,
    within_patient_variance = comp$gamma2 + comp$sigma2,
    total_sd = total_sd,
    lower_loa = mu - 1.96 * total_sd,
    upper_loa = mu + 1.96 * total_sd,
    log_start = opt$par,
    objective = opt$value,
    hessian_min_eigen = hessian_min_eigen,
    clustering_structure = "patient_and_stay_nested"
  )
}

reml_objective <- function(log_par, n, mean_y, within_ss, freq) {
  tau2 <- exp(log_par[[1L]])
  sigma2 <- exp(log_par[[2L]])
  denom <- sigma2 + n * tau2
  info <- sum(freq * n / denom)
  if (!is.finite(info) || info <= 0) return(.Machine$double.xmax / 100)
  mu <- sum(freq * n * mean_y / denom) / info
  logdet <- sum(freq * ((n - 1) * log(sigma2) + log(denom)))
  quad <- sum(freq * within_ss / sigma2) + sum(freq * n * (mean_y - mu)^2 / denom)
  n_total <- sum(freq * n)
  if (n_total <= 1) return(.Machine$double.xmax / 100)
  0.5 * (logdet + log(info) + quad + (n_total - 1) * log(2 * pi))
}

fit_reml_sufficient <- function(cs, freq = rep.int(1, nrow(cs)), start = NULL, require_hessian = FALSE) {
  if (length(freq) != nrow(cs)) stop("FREQUENCY_LENGTH_MISMATCH", call. = FALSE)
  keep <- freq > 0
  n <- cs$n[keep]
  mean_y <- cs$mean_y[keep]
  within_ss <- cs$within_ss[keep]
  freq <- as.numeric(freq[keep])
  if (sum(freq) < 2 || sum(freq * n) < 3) stop("INSUFFICIENT_CLUSTERS", call. = FALSE)
  if (is.null(start)) {
    total_mean <- sum(freq * n * mean_y) / sum(freq * n)
    total_ss <- sum(freq * within_ss) + sum(freq * n * (mean_y - total_mean)^2)
    total_var <- max(total_ss / max(sum(freq * n) - 1, 1), 1e-4)
    sigma2_0 <- max(sum(freq * within_ss) / max(sum(freq * (n - 1)), 1), total_var * 0.25, 1e-5)
    tau2_0 <- max(total_var - sigma2_0, total_var * 0.1, 1e-5)
    start <- log(c(tau2_0, sigma2_0))
  }
  opt <- optim(
    par = start,
    fn = reml_objective,
    n = n,
    mean_y = mean_y,
    within_ss = within_ss,
    freq = freq,
    method = "L-BFGS-B",
    lower = log(c(1e-10, 1e-10)),
    upper = log(c(1e4, 1e4)),
    control = list(maxit = 300, factr = 1e7)
  )
  if (opt$convergence != 0L || !is.finite(opt$value)) stop("REML_OPTIMIZATION_FAILED", call. = FALSE)
  tau2 <- exp(opt$par[[1L]])
  sigma2 <- exp(opt$par[[2L]])
  denom <- sigma2 + n * tau2
  info <- sum(freq * n / denom)
  mu <- sum(freq * n * mean_y / denom) / info
  total_sd <- sqrt(tau2 + sigma2)
  hessian_min_eigen <- NA_real_
  if (require_hessian) {
    h <- optimHess(opt$par, reml_objective, n = n, mean_y = mean_y, within_ss = within_ss, freq = freq)
    ev <- eigen((h + t(h)) / 2, symmetric = TRUE, only.values = TRUE)$values
    hessian_min_eigen <- min(ev)
    if (any(!is.finite(ev)) || hessian_min_eigen <= 0) stop("REML_HESSIAN_NOT_POSITIVE", call. = FALSE)
    if (tau2 <= 1e-8 || sigma2 <= 1e-8) stop("REML_SINGULAR_VARIANCE", call. = FALSE)
  }
  list(
    mean_bias = mu,
    between_patient_variance = tau2,
    between_stay_within_patient_variance = NA_real_,
    within_stay_residual_variance = sigma2,
    within_patient_variance = sigma2,
    total_sd = total_sd,
    lower_loa = mu - 1.96 * total_sd,
    upper_loa = mu + 1.96 * total_sd,
    log_start = opt$par,
    objective = opt$value,
    hessian_min_eigen = hessian_min_eigen,
    clustering_structure = "patient_random_intercept"
  )
}

weighted_quantiles_from_counts <- function(values, counts, probs = c(0.025, 0.975)) {
  total <- sum(counts)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(probs)))
  cumulative <- cumsum(counts) / total
  vapply(probs, function(p) values[[which(cumulative >= p)[1L]]], numeric(1L))
}

build_histogram_matrix <- function(z, cs) {
  values <- sort(unique(z$bias))
  h <- sparseMatrix(
    i = z$patient_index,
    j = match(z$bias, values),
    x = 1,
    dims = c(nrow(cs), length(values))
  )
  list(values = values, counts = h)
}

original_distribution_metrics <- function(z, cs, hist) {
  empirical <- weighted_quantiles_from_counts(hist$values, as.numeric(colSums(hist$counts)))
  mean_pair <- mean(z$bias)
  arms_pair <- sqrt(mean(z$bias^2))
  patient_mean <- mean(cs$mean_y)
  patient_second <- mean(cs$mean_y2)
  patient_sd <- sqrt(max(patient_second - patient_mean^2, 0))
  patient_counts <- as.numeric(crossprod(1 / cs$n, hist$counts))
  patient_empirical <- weighted_quantiles_from_counts(hist$values, patient_counts)
  list(
    pair_mean = mean_pair,
    pair_arms = arms_pair,
    empirical_p2_5 = empirical[[1L]],
    empirical_p97_5 = empirical[[2L]],
    balanced_mean = patient_mean,
    balanced_sd = patient_sd,
    balanced_lower_loa = patient_mean - 1.96 * patient_sd,
    balanced_upper_loa = patient_mean + 1.96 * patient_sd,
    balanced_arms = sqrt(patient_second),
    balanced_empirical_p2_5 = patient_empirical[[1L]],
    balanced_empirical_p97_5 = patient_empirical[[2L]]
  )
}

cluster_bootstrap <- function(z, cs, fit, b = 2000L, seed = 1L, include_histograms = TRUE) {
  set.seed(seed)
  m <- nrow(cs)
  hist <- if (include_histograms) build_histogram_matrix(z, cs) else NULL
  out <- matrix(NA_real_, nrow = b, ncol = 14L)
  colnames(out) <- c(
    "mean_bias", "lower_loa", "upper_loa", "between_variance", "within_variance",
    "pair_arms", "empirical_p2_5", "empirical_p97_5",
    "balanced_mean", "balanced_lower_loa", "balanced_upper_loa", "balanced_arms",
    "balanced_empirical_p2_5", "balanced_empirical_p97_5"
  )
  failure_reason <- character(b)
  for (r in seq_len(b)) {
    freq <- tabulate(sample.int(m, m, replace = TRUE), nbins = m)
    current <- try(fit_reml_sufficient(cs, freq = freq, start = fit$log_start, require_hessian = FALSE), silent = TRUE)
    if (inherits(current, "try-error")) {
      failure_reason[[r]] <- conditionMessage(attr(current, "condition"))
      next
    }
    pair_n <- sum(freq * cs$n)
    pair_second <- sum(freq * cs$n * cs$mean_y2) / pair_n
    out[r, "mean_bias"] <- current$mean_bias
    out[r, "lower_loa"] <- current$lower_loa
    out[r, "upper_loa"] <- current$upper_loa
    out[r, "between_variance"] <- current$between_patient_variance
    out[r, "within_variance"] <- current$within_patient_variance
    out[r, "pair_arms"] <- sqrt(pair_second)
    patient_mean <- sum(freq * cs$mean_y) / m
    patient_second <- sum(freq * cs$mean_y2) / m
    patient_sd <- sqrt(max(patient_second - patient_mean^2, 0))
    out[r, "balanced_mean"] <- patient_mean
    out[r, "balanced_lower_loa"] <- patient_mean - 1.96 * patient_sd
    out[r, "balanced_upper_loa"] <- patient_mean + 1.96 * patient_sd
    out[r, "balanced_arms"] <- sqrt(patient_second)
    if (include_histograms) {
      pair_counts <- as.numeric(crossprod(freq, hist$counts))
      q <- weighted_quantiles_from_counts(hist$values, pair_counts)
      out[r, "empirical_p2_5"] <- q[[1L]]
      out[r, "empirical_p97_5"] <- q[[2L]]
      balanced_counts <- as.numeric(crossprod(freq / cs$n, hist$counts))
      bq <- weighted_quantiles_from_counts(hist$values, balanced_counts)
      out[r, "balanced_empirical_p2_5"] <- bq[[1L]]
      out[r, "balanced_empirical_p97_5"] <- bq[[2L]]
    }
  }
  success <- complete.cases(out[, c("mean_bias", "lower_loa", "upper_loa")])
  list(values = out, success = success, failures = failure_reason, seed = seed)
}

nested_bootstrap_worker <- function(r) {
  set.seed(.stage403_replicate_seeds[[r]])
  freq <- tabulate(sample.int(.stage403_m, .stage403_m, replace = TRUE), nbins = .stage403_m)
  current <- try(fit_nested_reml_sufficient(
    .stage403_ss, patient_count = .stage403_m, freq = freq,
    start = .stage403_fit_start, require_hessian = FALSE
  ), silent = TRUE)
  if (inherits(current, "try-error")) {
    return(list(values = rep(NA_real_, 14L), failure = conditionMessage(attr(current, "condition"))))
  }
  values <- rep(NA_real_, 14L)
  names(values) <- .stage403_output_names
  pair_n <- sum(freq * .stage403_cs$n)
  pair_second <- sum(freq * .stage403_cs$n * .stage403_cs$mean_y2) / pair_n
  values["mean_bias"] <- current$mean_bias
  values["lower_loa"] <- current$lower_loa
  values["upper_loa"] <- current$upper_loa
  values["between_variance"] <- current$between_patient_variance
  values["within_variance"] <- current$within_patient_variance
  values["pair_arms"] <- sqrt(pair_second)
  patient_mean <- sum(freq * .stage403_cs$mean_y) / .stage403_m
  patient_second <- sum(freq * .stage403_cs$mean_y2) / .stage403_m
  patient_sd <- sqrt(max(patient_second - patient_mean^2, 0))
  values["balanced_mean"] <- patient_mean
  values["balanced_lower_loa"] <- patient_mean - 1.96 * patient_sd
  values["balanced_upper_loa"] <- patient_mean + 1.96 * patient_sd
  values["balanced_arms"] <- sqrt(patient_second)
  if (.stage403_include_histograms) {
    pair_counts <- as.numeric(crossprod(freq, .stage403_hist$counts))
    q <- weighted_quantiles_from_counts(.stage403_hist$values, pair_counts)
    values["empirical_p2_5"] <- q[[1L]]
    values["empirical_p97_5"] <- q[[2L]]
    balanced_counts <- as.numeric(crossprod(freq / .stage403_cs$n, .stage403_hist$counts))
    bq <- weighted_quantiles_from_counts(.stage403_hist$values, balanced_counts)
    values["balanced_empirical_p2_5"] <- bq[[1L]]
    values["balanced_empirical_p97_5"] <- bq[[2L]]
  }
  list(values = values, failure = "")
}

cluster_bootstrap_nested <- function(z, cs, ss, fit, b = 2000L, seed = 1L, include_histograms = TRUE) {
  output_names <- c(
    "mean_bias", "lower_loa", "upper_loa", "between_variance", "within_variance",
    "pair_arms", "empirical_p2_5", "empirical_p97_5",
    "balanced_mean", "balanced_lower_loa", "balanced_upper_loa", "balanced_arms",
    "balanced_empirical_p2_5", "balanced_empirical_p97_5"
  )
  set.seed(seed)
  replicate_seeds <- sample.int(.Machine$integer.max, b, replace = FALSE)
  m <- nrow(cs)
  hist <- if (include_histograms) build_histogram_matrix(z, cs) else NULL
  worker_environment <- new.env(parent = environment())
  worker_environment$.stage403_replicate_seeds <- replicate_seeds
  worker_environment$.stage403_m <- m
  worker_environment$.stage403_ss <- ss
  worker_environment$.stage403_cs <- cs
  worker_environment$.stage403_fit_start <- fit$log_start
  worker_environment$.stage403_output_names <- output_names
  worker_environment$.stage403_include_histograms <- include_histograms
  worker_environment$.stage403_hist <- hist
  worker_environment$nested_bootstrap_worker <- nested_bootstrap_worker
  worker_environment$fit_nested_reml_sufficient <- fit_nested_reml_sufficient
  worker_environment$nested_reml_objective <- nested_reml_objective
  worker_environment$nested_reml_components <- nested_reml_components
  worker_environment$weighted_quantiles_from_counts <- weighted_quantiles_from_counts
  cores <- if (b >= 100L) min(4L, max(1L, parallel::detectCores(logical = FALSE)), b) else 1L
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    export_names <- ls(worker_environment, all.names = TRUE)
    parallel::clusterExport(cl, export_names, envir = worker_environment)
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)); suppressPackageStartupMessages(library(Matrix)); NULL })
    result <- parallel::parLapplyLB(cl, seq_len(b), nested_bootstrap_worker)
  } else {
    local_worker <- nested_bootstrap_worker
    environment(local_worker) <- worker_environment
    result <- lapply(seq_len(b), local_worker)
  }
  out <- do.call(rbind, lapply(result, `[[`, "values"))
  colnames(out) <- output_names
  failure_reason <- vapply(result, `[[`, character(1L), "failure")
  success <- complete.cases(out[, c("mean_bias", "lower_loa", "upper_loa")])
  list(values = out, success = success, failures = failure_reason, seed = seed, parallel_cores = cores)
}

bootstrap_interval <- function(boot, variable) {
  x <- boot$values[boot$success, variable]
  if (!length(x)) return(c(NA_real_, NA_real_))
  as.numeric(quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 7))
}

summarize_overall <- function(z, b = 2000L, seed = 1L, nested = FALSE) {
  cs <- cluster_summaries(copy(z))
  ss <- if (nested) nested_stay_summaries(copy(z)) else NULL
  fit <- if (nested) fit_nested_reml_sufficient(ss, nrow(cs), require_hessian = TRUE) else fit_reml_sufficient(cs, require_hessian = TRUE)
  hist <- build_histogram_matrix(z, cs)
  dist <- original_distribution_metrics(z, cs, hist)
  boot <- if (nested) cluster_bootstrap_nested(z, cs, ss, fit, b = b, seed = seed, include_histograms = TRUE) else cluster_bootstrap(z, cs, fit, b = b, seed = seed, include_histograms = TRUE)
  if (sum(boot$success) < ceiling(0.95 * b)) stop("BOOTSTRAP_SUCCESS_BELOW_95_PERCENT", call. = FALSE)
  ci <- lapply(colnames(boot$values), function(v) bootstrap_interval(boot, v))
  names(ci) <- colnames(boot$values)
  list(fit = fit, dist = dist, boot = boot, ci = ci, clusters = cs, stays = ss)
}

summarize_stratum <- function(z, b = 2000L, seed = 1L, nested = FALSE) {
  z <- copy(z)
  z[, patient_index := match(patient, unique(patient))]
  cs <- cluster_summaries(z)
  ss <- if (nested) nested_stay_summaries(copy(z)) else NULL
  fit <- if (nested) fit_nested_reml_sufficient(ss, nrow(cs), require_hessian = FALSE) else fit_reml_sufficient(cs, require_hessian = FALSE)
  boot <- if (nested) cluster_bootstrap_nested(z, cs, ss, fit, b = b, seed = seed, include_histograms = FALSE) else cluster_bootstrap(z, cs, fit, b = b, seed = seed, include_histograms = FALSE)
  if (sum(boot$success) < ceiling(0.95 * b)) stop("STRATUM_BOOTSTRAP_SUCCESS_BELOW_95_PERCENT", call. = FALSE)
  list(fit = fit, boot = boot, ci = lapply(c("mean_bias", "lower_loa", "upper_loa"), function(v) bootstrap_interval(boot, v)))
}

fit_proportional_bias <- function(z, nested = FALSE) {
  suppressPackageStartupMessages(library(lme4))
  formula <- if (nested) bias ~ centred_mean + (1 | patient) + (1 | stay) else bias ~ centred_mean + (1 | patient)
  model_control <- if (nested) {
    lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 500000), calc.derivs = TRUE)
  } else {
    lmerControl(optimizer = "nloptwrap", calc.derivs = TRUE)
  }
  model <- lmer(
    formula, data = z, REML = TRUE,
    control = model_control
  )
  singular_fit <- isSingular(model, tol = 1e-8)
  opt <- model@optinfo
  if (!is.null(opt$conv$lme4$messages)) stop("PROPORTIONAL_MODEL_CONVERGENCE_MESSAGE", call. = FALSE)
  co <- summary(model)$coefficients
  intercept_at_90 <- co["(Intercept)", "Estimate"]
  slope <- co["centred_mean", "Estimate"]
  se <- co["centred_mean", "Std. Error"]
  list(model = model, intercept_at_90 = intercept_at_90, slope = slope, se = se,
       lower = slope - 1.96 * se, upper = slope + 1.96 * se,
       singular_fit = singular_fit)
}

fit_heteroscedasticity <- function(z, nested = FALSE) {
  suppressPackageStartupMessages(library(nlme))
  ctl <- lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50, opt = "nlminb", returnObject = FALSE)
  random_formula <- if (nested) ~1 | patient/stay else ~1 | patient
  constant <- lme(
    bias ~ centred_mean, random = random_formula, data = z,
    method = "ML", control = ctl, na.action = na.fail
  )
  varying <- lme(
    bias ~ centred_mean, random = random_formula,
    weights = varExp(form = ~ centred_mean), data = z,
    method = "ML", control = ctl, na.action = na.fail
  )
  comparison <- anova(constant, varying)
  delta <- as.numeric(coef(varying$modelStruct$varStruct, unconstrained = FALSE))
  interval_try <- try(intervals(varying, which = "var-cov")$varStruct, silent = TRUE)
  ci_method <- "approximate_variance_covariance"
  if (!inherits(interval_try, "try-error")) {
    delta_lower <- as.numeric(interval_try[1L, "lower"])
    delta_upper <- as.numeric(interval_try[1L, "upper"])
  } else {
    ci_method <- "profile_likelihood_fallback"
    ll_max <- as.numeric(logLik(varying))
    cutoff <- qchisq(0.95, df = 1)
    profile_objective <- function(candidate_delta) {
      fixed_fit <- try(lme(
        bias ~ centred_mean, random = random_formula,
        weights = varExp(form = ~ centred_mean, fixed = candidate_delta), data = z,
        method = "ML", control = ctl, na.action = na.fail
      ), silent = TRUE)
      if (inherits(fixed_fit, "try-error")) return(NA_real_)
      2 * (ll_max - as.numeric(logLik(fixed_fit))) - cutoff
    }
    find_endpoint <- function(direction) {
      inner <- delta
      f_inner <- -cutoff
      step <- 0.01
      for (iteration in seq_len(12L)) {
        outer <- delta + direction * step
        f_outer <- profile_objective(outer)
        if (is.finite(f_outer) && f_outer >= 0) {
          interval <- sort(c(inner, outer))
          return(uniroot(function(d) profile_objective(d), interval = interval, tol = 1e-6)$root)
        }
        if (is.finite(f_outer)) {
          inner <- outer
          f_inner <- f_outer
        }
        step <- step * 1.8
      }
      stop("HETEROSCEDASTIC_PROFILE_CI_ENDPOINT_FAILED", call. = FALSE)
    }
    delta_lower <- find_endpoint(-1)
    delta_upper <- find_endpoint(1)
  }
  list(
    constant = constant,
    varying = varying,
    delta = delta,
    delta_lower = delta_lower,
    delta_upper = delta_upper,
    ci_method = ci_method,
    likelihood_ratio = as.numeric(comparison$L.Ratio[[2L]]),
    p_value = as.numeric(comparison$`p-value`[[2L]]),
    aic_constant = AIC(constant),
    aic_varying = AIC(varying),
    residual_sigma_at_90 = sigma(varying)
  )
}

descriptive_smooth <- function(z, span = 0.75, grid = NULL) {
  if (is.null(grid)) {
    lower <- max(70, floor(min(z$paired_mean, na.rm = TRUE) * 10) / 10)
    upper <- min(100, ceiling(max(z$paired_mean, na.rm = TRUE) * 10) / 10)
    grid <- seq(lower, upper, by = 0.1)
  }
  fit <- loess(
    bias ~ paired_mean, data = z, span = span, degree = 2,
    family = "gaussian", surface = "direct",
    control = loess.control(trace.hat = "approximate")
  )
  out <- data.table(paired_mean = grid, smooth_bias = as.numeric(predict(fit, newdata = data.frame(paired_mean = grid))))
  out[is.finite(smooth_bias)]
}
