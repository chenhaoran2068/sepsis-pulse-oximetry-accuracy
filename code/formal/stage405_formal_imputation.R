s405f_run_imputation <- function(x, cohort, eligible_fields, m = 30L, maxit = 20L,
                                 seed = 20260907L, print_flag = FALSE) {
  s405m_stop_if(m != 30L || maxit != 20L, "S405F_FORMAL_DESIGN_MUST_BE_30X20")
  prepared <- s405m_prepare_analytic_states(x, cohort)
  encoded <- s405m_encode_imputation_data(prepared$data, cohort, eligible_fields)
  stay_result <- s405m_run_stay_imputation(
    encoded, cohort, as.integer(m), as.integer(maxit), as.integer(seed), print_flag
  )
  raw <- data.table::as.data.table(encoded$data)
  pair_results <- vector("list", as.integer(m))
  completed <- vector("list", as.integer(m))
  for (i in seq_len(m)) {
    stay_completed <- data.table::as.data.table(stay_result$completed[[i]])
    stay_values <- stay_completed[, c("cluster_patient", "cluster_stay", encoded$metadata$stay_fields), with = FALSE]
    s405m_assert_unique(stay_values, c("cluster_patient", "cluster_stay"), "S405F_COMPLETED_STAY_COMPOSITE_KEY_DUPLICATE")
    pair_data <- merge(
      raw[, setdiff(names(raw), encoded$metadata$stay_fields), with = FALSE],
      stay_values, by = c("cluster_patient", "cluster_stay"), all.x = TRUE, sort = FALSE
    )
    s405m_stop_if(nrow(pair_data) != nrow(raw), "S405F_STAY_MAPPING_CARDINALITY")
    pair_data <- as.data.frame(pair_data[, names(encoded$data), with = FALSE])
    pair_results[[i]] <- s405m_run_pair_chain(
      pair_data, cohort, encoded$metadata$pair_fields,
      maxit = as.integer(maxit), seed = as.integer(seed + 10000L + i),
      print_flag = print_flag, chain = i
    )
    completed[[i]] <- s405m_complete_data(
      NULL,
      list(data = pair_results[[i]]$completed, metadata = encoded$metadata),
      1L
    )
  }
  logged <- data.table::rbindlist(c(
    list(stay_result$logged_events),
    lapply(pair_results, function(result) result$logged_events)
  ), fill = TRUE)
  trace <- data.table::rbindlist(c(
    list(stay_result$trace),
    lapply(pair_results, function(result) result$trace)
  ), fill = TRUE)
  predictor_ledger <- data.table::rbindlist(c(
    list(stay_result$predictor_ledger),
    lapply(pair_results, function(result) result$predictor_ledger)
  ), fill = TRUE)
  sex_summary <- prepared$sex_ledger[, .(
    stay_n = .N,
    source_observed_n = sum(source_observed),
    patient_reconciled_n = sum(reconciled_from_patient),
    analytic_observed_before_imputation_n = sum(analytic_observed),
    residual_missing_n = sum(residual_missing)
  )]
  convergence <- s405c_run_adapter(trace, 30L, 20L, "formal")
  list(
    stay_mids = stay_result$mids,
    pair_mids = lapply(pair_results, function(result) result$mids),
    completed = completed,
    encoded = encoded,
    spec = list(
      stay_targets = stay_result$targets,
      pair_targets = unique(unlist(lapply(pair_results, function(result) result$targets))),
      pair_imputation_structure = if (cohort == "eICU") "pairs_within_stay_with_patient_shared_auxiliaries" else "pairs_within_patient_stay",
      patient_auxiliary_n = max(vapply(pair_results, function(result) result$patient_auxiliary_n, integer(1L))),
      patient_auxiliary_source_fields = unique(unlist(lapply(pair_results, function(result) result$patient_auxiliary_source_fields))),
      pair_nonimputed_missing_predictors_excluded = unique(unlist(lapply(pair_results, function(result) result$nonimputed_missing_predictors_excluded))),
      pair_redundant_stay_auxiliaries_excluded = unique(unlist(lapply(pair_results, function(result) result$redundant_stay_auxiliaries_excluded))),
      formal_convergence_adapter_applied = TRUE,
      formal_convergence_screen_pass = convergence$gate$formal_convergence_screen_pass[[1L]],
      convergence_interpretation = convergence$gate$interpretation[[1L]]
    ),
    sex_reconciliation_summary = sex_summary,
    predictor_ledger = predictor_ledger,
    logged_events = logged,
    trace = trace,
    convergence = convergence,
    m = as.integer(m),
    maxit = as.integer(maxit)
  )
}
