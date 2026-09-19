s405c_stop_if <- function(condition, code) {
  if (isTRUE(condition)) stop(code, call. = FALSE)
  invisible(TRUE)
}

s405c_required_trace_fields <- function() {
  c("imputation_block", "feature", "iteration", "chain", "chain_mean", "chain_variance")
}

s405c_validate_trace <- function(trace, expected_chains, expected_iterations) {
  trace <- data.table::as.data.table(trace)
  required <- s405c_required_trace_fields()
  s405c_stop_if(!all(required %in% names(trace)), "S405C_TRACE_SCHEMA_MISSING")
  s405c_stop_if(!nrow(trace), "S405C_TRACE_EMPTY")
  s405c_stop_if(anyNA(trace[, c("imputation_block", "feature", "iteration", "chain", "chain_mean"), with = FALSE]), "S405C_TRACE_HAS_MISSING_REQUIRED_VALUE")
  s405c_stop_if(any(!is.finite(trace$chain_mean)), "S405C_CHAIN_MEAN_NONFINITE")
  variance_state <- trace[, .(
    variance_missing_n = sum(is.na(chain_variance)),
    variance_nonmissing_n = sum(!is.na(chain_variance)),
    variance_nonfinite_n = sum(!is.na(chain_variance) & !is.finite(chain_variance)),
    variance_negative_n = sum(!is.na(chain_variance) & chain_variance < 0)
  ), by = .(imputation_block, feature)]
  s405c_stop_if(any(variance_state$variance_missing_n > 0L & variance_state$variance_nonmissing_n > 0L), "S405C_CHAIN_VARIANCE_PARTIALLY_MISSING")
  s405c_stop_if(any(variance_state$variance_nonfinite_n > 0L), "S405C_CHAIN_VARIANCE_NONFINITE")
  s405c_stop_if(any(variance_state$variance_negative_n > 0L), "S405C_CHAIN_VARIANCE_NEGATIVE")
  duplicate_n <- trace[, .N, by = .(imputation_block, feature, iteration, chain)][N != 1L, .N]
  s405c_stop_if(duplicate_n > 0L, "S405C_TRACE_DUPLICATE_CELL")
  completeness <- trace[, .(
    chain_n = data.table::uniqueN(chain),
    iteration_n = data.table::uniqueN(iteration),
    cell_n = .N,
    chain_min = min(chain), chain_max = max(chain),
    iteration_min = min(iteration), iteration_max = max(iteration)
  ), by = .(imputation_block, feature)]
  complete <- completeness[
    chain_n == expected_chains & iteration_n == expected_iterations &
      cell_n == expected_chains * expected_iterations &
      chain_min == 1L & chain_max == expected_chains &
      iteration_min == 1L & iteration_max == expected_iterations
  ]
  s405c_stop_if(nrow(complete) != nrow(completeness), "S405C_TRACE_CARDINALITY_MISMATCH")
  completeness
}

s405c_trace_matrix <- function(trace, block, feature, parameter) {
  feature_value <- feature
  x <- data.table::copy(trace[imputation_block == block & get("feature") == feature_value])
  x[, value := if (parameter == "mean") chain_mean else sqrt(chain_variance)]
  wide <- data.table::dcast(x, iteration ~ chain, value.var = "value")
  data.table::setorder(wide, iteration)
  as.matrix(wide[, -"iteration"])
}

s405c_safe_rhat <- function(matrix_value) {
  if (length(unique(as.numeric(matrix_value))) <= 1L) return(1)
  value <- suppressWarnings(posterior::rhat(matrix_value))
  if (!is.finite(value)) NA_real_ else as.numeric(value)
}

s405c_lag1 <- function(matrix_value) {
  if (nrow(matrix_value) < 3L) return(NA_real_)
  ac <- vapply(seq_len(ncol(matrix_value)), function(chain) {
    value <- matrix_value[, chain]
    if (stats::sd(value) == 0) return(0)
    suppressWarnings(stats::cor(value[-length(value)], value[-1L]))
  }, numeric(1L))
  mean(abs(ac), na.rm = TRUE)
}

s405c_half_drift <- function(matrix_value) {
  if (length(unique(as.numeric(matrix_value))) <= 1L) return(0)
  split <- floor(nrow(matrix_value) / 2L)
  early <- matrix_value[seq_len(split), , drop = FALSE]
  late <- matrix_value[seq.int(split + 1L, nrow(matrix_value)), , drop = FALSE]
  denominator <- stats::sd(as.numeric(matrix_value))
  if (!is.finite(denominator) || denominator == 0) return(0)
  (mean(late) - mean(early)) / denominator
}

s405c_parameter_diagnostic <- function(trace, block, feature, parameter) {
  matrix_value <- s405c_trace_matrix(trace, block, feature, parameter)
  applicable <- !all(is.na(matrix_value))
  if (!applicable) {
    summary <- data.table::data.table(
      imputation_block = block, feature = feature, parameter = parameter,
      parameter_applicable = FALSE,
      chain_n = ncol(matrix_value), iteration_n = nrow(matrix_value),
      constant_trajectory = NA,
      rank_normalized_split_rhat = NA_real_,
      mean_absolute_lag1_autocorrelation = NA_real_,
      standardized_half_drift = NA_real_,
      review_required = FALSE,
      blocking_failure = FALSE,
      blocking_reason = "",
      not_applicable_reason = "chain_variance_undefined_insufficient_imputed_cells"
    )
    trajectory <- data.table::data.table(
      imputation_block = block, feature = feature, parameter = parameter,
      iteration = seq_len(nrow(matrix_value)),
      across_chain_min = NA_real_, across_chain_median = NA_real_,
      across_chain_max = NA_real_, cumulative_rhat = NA_real_,
      transition_autocorrelation = NA_real_
    )
    return(list(summary = summary, trajectory = trajectory))
  }
  constant <- length(unique(as.numeric(matrix_value))) <= 1L
  rhat <- s405c_safe_rhat(matrix_value)
  lag1 <- s405c_lag1(matrix_value)
  drift <- s405c_half_drift(matrix_value)
  blocking <- !is.finite(rhat) || !is.finite(lag1) || !is.finite(drift) ||
    rhat > 1.10 || lag1 > 0.80 || abs(drift) > 0.50
  review <- blocking || rhat > 1.05 || lag1 > 0.50 || abs(drift) > 0.25
  summary <- data.table::data.table(
    imputation_block = block, feature = feature, parameter = parameter,
    parameter_applicable = TRUE,
    chain_n = ncol(matrix_value), iteration_n = nrow(matrix_value),
    constant_trajectory = constant,
    rank_normalized_split_rhat = rhat,
    mean_absolute_lag1_autocorrelation = lag1,
    standardized_half_drift = drift,
    review_required = review,
    blocking_failure = blocking,
    blocking_reason = paste(c(
      if (!is.finite(rhat)) "rhat_nonfinite" else if (rhat > 1.10) "rhat_gt_1.10" else NULL,
      if (!is.finite(lag1)) "lag1_nonfinite" else if (lag1 > 0.80) "lag1_gt_0.80" else NULL,
      if (!is.finite(drift)) "drift_nonfinite" else if (abs(drift) > 0.50) "absolute_drift_gt_0.50" else NULL
    ), collapse = "|"),
    not_applicable_reason = ""
  )
  cumulative <- data.table::rbindlist(lapply(seq_len(nrow(matrix_value)), function(iteration) {
    current <- matrix_value[seq_len(iteration), , drop = FALSE]
    transition_ac <- if (iteration < 2L) NA_real_ else {
      a <- current[iteration - 1L, ]
      b <- current[iteration, ]
      if (stats::sd(a) == 0 || stats::sd(b) == 0) 0 else suppressWarnings(stats::cor(a, b))
    }
    data.table::data.table(
      imputation_block = block, feature = feature, parameter = parameter,
      iteration = iteration,
      across_chain_min = min(current[iteration, ]),
      across_chain_median = stats::median(current[iteration, ]),
      across_chain_max = max(current[iteration, ]),
      cumulative_rhat = if (iteration < 4L) NA_real_ else s405c_safe_rhat(current),
      transition_autocorrelation = transition_ac
    )
  }))
  list(summary = summary, trajectory = cumulative)
}

s405c_run_adapter <- function(trace, expected_chains, expected_iterations,
                              mode = c("pilot_interface", "formal")) {
  mode <- match.arg(mode)
  trace <- data.table::as.data.table(trace)
  completeness <- s405c_validate_trace(trace, expected_chains, expected_iterations)
  formal_design <- expected_chains == 30L && expected_iterations == 20L
  s405c_stop_if(mode == "formal" && !formal_design, "S405C_FORMAL_DESIGN_MUST_BE_30X20")
  keys <- unique(trace[, .(imputation_block, feature)])
  diagnostics <- lapply(seq_len(nrow(keys)), function(i) {
    lapply(c("mean", "sd"), function(parameter) {
      s405c_parameter_diagnostic(trace, keys$imputation_block[[i]], keys$feature[[i]], parameter)
    })
  })
  flat <- unlist(diagnostics, recursive = FALSE)
  summary <- data.table::rbindlist(lapply(flat, `[[`, "summary"), fill = TRUE)
  trajectory <- data.table::rbindlist(lapply(flat, `[[`, "trajectory"), fill = TRUE)
  formal_gate_eligible <- mode == "formal" && formal_design
  list(
    completeness = completeness,
    summary = summary,
    trajectory = trajectory,
    gate = data.table::data.table(
      mode = mode,
      expected_chains = expected_chains,
      expected_iterations = expected_iterations,
      formal_gate_eligible = formal_gate_eligible,
      diagnostic_row_n = nrow(summary),
      review_required_n = sum(summary$review_required),
      blocking_failure_n = sum(summary$blocking_failure),
      formal_convergence_screen_pass = formal_gate_eligible && !any(summary$blocking_failure),
      interpretation = if (formal_gate_eligible) {
        "formal_computational_screen_requires_trajectory_review"
      } else {
        "interface_validation_only_not_formal_convergence_evidence"
      }
    )
  )
}
