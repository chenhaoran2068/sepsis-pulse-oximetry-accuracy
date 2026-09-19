args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: ANALYSIS_INPUT_DIR MODEL_OUTPUT_DIR NEW_QA_DIR", call. = FALSE)
input <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
model <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
qa <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
if (file.exists(qa) || dir.exists(qa)) stop("QA_OUTPUT_EXISTS_REFUSE_OVERWRITE", call. = FALSE)
suppressPackageStartupMessages(library(digest))
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
checks <- list()
check <- function(name, condition) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(condition))
}
aggregate <- file.path(model, "aggregate")
summary <- read.delim(file.path(aggregate, "run_summary.tsv"), check.names = FALSE)
manifest <- read.delim(file.path(aggregate, "completed_dataset_manifest.tsv"), check.names = FALSE)
gate <- read.delim(file.path(aggregate, "convergence_gate.tsv"), check.names = FALSE)
variable <- read.delim(file.path(aggregate, "imputation_variable_validation.tsv"), check.names = FALSE)
constancy <- read.delim(file.path(aggregate, "stay_value_constancy.tsv"), check.names = FALSE)
sex <- read.delim(file.path(aggregate, "patient_sex_constancy.tsv"), check.names = FALSE)
domain <- read.delim(file.path(aggregate, "imputed_value_domain_audit.tsv"), check.names = FALSE)
primary <- read.delim(file.path(aggregate, "lmm_primary_pooled.tsv"), check.names = FALSE)
balanced <- read.delim(file.path(aggregate, "lmm_patient_balanced_pooled.tsv"), check.names = FALSE)
check("candidate_status", identical(readLines(file.path(model, "status.txt"), warn = FALSE),
                                    "FORMAL_MI_LMM_CANDIDATE_COMPLETE_QA_PENDING"))
check("one_cohort_window", nrow(summary) == 1L && summary$window_minutes[[1L]] %in% c(60L, 5L))
source <- readRDS(file.path(input, paste0("lmm_input_", summary$window_minutes[[1L]], ".rds")))
check("input_scope", all(source$cohort == summary$cohort[[1L]]) &&
                       all(source$pair_window_minutes == summary$window_minutes[[1L]]))
check("input_hash", identical(as.character(summary$input_sha256[[1L]]),
          digest::digest(file.path(input, paste0("lmm_input_", summary$window_minutes[[1L]], ".rds")),
                         file = TRUE, algo = "sha256")))
check("count", nrow(source) == summary$pair_n[[1L]] &&
               length(unique(source$patient_key_internal)) == summary$patient_n[[1L]])
check("exact_30x20", summary$imputations[[1L]] == 30L && summary$iterations[[1L]] == 20L &&
       nrow(manifest) == 30L && identical(sort(as.integer(manifest$imputation)), seq_len(30L)))
check("convergence_gate", nrow(gate) == 1L && gate$formal_gate_eligible[[1L]] &&
       gate$expected_chains[[1L]] == 30L && gate$expected_iterations[[1L]] == 20L &&
       gate$blocking_failure_n[[1L]] == 0L)
check("variable_validation", nrow(variable) > 0L && all(variable$observed_unchanged) &&
       all(variable$completed_no_missing))
check("stay_constancy", nrow(constancy) > 0L && all(constancy$stay_constant))
check("patient_sex_constancy", nrow(sex) == 0L || all(sex$patient_constant))
check("imputed_domain", nrow(domain) > 0L && all(domain$domain_pass))
check("no_individual_saturation_terms", !any(grepl("^spo2|^sao2", primary$term)) &&
       !any(grepl("^spo2|^sao2", balanced$term)))
check("oxygenation_background", "paired_mean_saturation_c90" %in% primary$term)
check("term_sets_match", setequal(primary$term, balanced$term))
check("finite_pooled", all(is.finite(primary$estimate)) && all(is.finite(primary$std_error)) &&
       all(is.finite(balanced$estimate)) && all(is.finite(balanced$std_error)))
check("ci_contains_estimate", all(primary$ci_lower <= primary$estimate & primary$estimate <= primary$ci_upper) &&
       all(balanced$ci_lower <= balanced$estimate & balanced$estimate <= balanced$ci_upper))
hash_ok <- vapply(seq_len(nrow(manifest)), function(i) {
  path <- file.path(model, "imputations", manifest$file[[i]])
  file.exists(path) && identical(as.character(manifest$sha256[[i]]),
                                 digest::digest(path, file = TRUE, algo = "sha256"))
}, logical(1))
check("all_imputation_hashes", all(hash_ok))
check("all_fit_checkpoints", all(file.exists(file.path(model, "checkpoints",
                              sprintf("fit_%02d.rds", seq_len(30L))))))
formal_gate <- read.delim(file.path(aggregate, "formal_lmm_gate.tsv"), check.names = FALSE)
fits <- lapply(file.path(model, "checkpoints", sprintf("fit_%02d.rds", seq_len(30L))), readRDS)
status_rows <- do.call(rbind, lapply(fits, function(x) x$primary$status))
design_rows <- do.call(rbind, lapply(fits, function(x) x$primary$design))
gvif_rows <- do.call(rbind, lapply(fits, function(x) x$primary$gvif))
residual_rows <- do.call(rbind, lapply(fits, function(x) x$primary$residual))
weight_rows <- do.call(rbind, lapply(fits, function(x) x$balanced$weight_status))
check("formal_gate_all_pass", nrow(formal_gate) == 11L && all(formal_gate$pass))
check("independent_lmm_convergence", nrow(status_rows) == 30L &&
      all(status_rows$converged & status_rows$hessian_positive_definite &
          status_rows$hessian_available & status_rows$opt_code == 0L))
check("independent_lmm_non_singular", !any(status_rows$singular) &&
      all(status_rows$random_variance_floor_pass))
check("independent_lmm_design_and_vif", all(design_rows$full_rank) &&
      all(gvif_rows$adjusted_gvif < 5))
check("independent_lmm_residuals", all(residual_rows$all_residuals_finite &
      residual_rows$all_fitted_finite))
check("independent_patient_balance", all(weight_rows$all_patient_total_weights_equal_one))
check("completed_row_counts", all(vapply(seq_len(nrow(manifest)), function(i) {
  path <- file.path(model, "imputations", manifest$file[[i]])
  nrow(readRDS(path)) == nrow(source) && manifest$row_n[[i]] == nrow(source)
}, logical(1))))
for (kind in c("primary", "balanced")) {
  table <- if (kind == "primary") primary else balanced
  rubin_ok <- vapply(seq_len(nrow(table)), function(j) {
    term <- table$term[[j]]
    coefficient <- variance <- numeric(30L)
    for (i in seq_len(30L)) {
      checkpoint <- readRDS(file.path(model, "checkpoints", sprintf("fit_%02d.rds", i)))[[kind]]
      index <- match(term, checkpoint$terms)
      if (is.na(index)) return(FALSE)
      coefficient[[i]] <- as.numeric(checkpoint$coefficients[[index]])
      variance[[i]] <- as.numeric(checkpoint$variances[[index]])
    }
    multiplier <- table$reporting_multiplier[[j]]
    estimate <- mean(coefficient) * multiplier
    standard_error <- sqrt(mean(variance) + (1 + 1 / 30) * stats::var(coefficient)) * abs(multiplier)
    isTRUE(all.equal(estimate, table$estimate[[j]], tolerance = 1e-8)) &&
      isTRUE(all.equal(standard_error, table$std_error[[j]], tolerance = 1e-8))
  }, logical(1))
  check(paste0("independent_rubin_", kind), all(rubin_ok))
}
results <- do.call(rbind, checks)
utils::write.table(results, file.path(qa, "checks.tsv"), sep = "\t", row.names = FALSE,
                   quote = FALSE)
status <- if (all(results$pass)) "FORMAL_MI_LMM_INDEPENDENT_QA_PASS" else
  "FORMAL_MI_LMM_INDEPENDENT_QA_FAIL"
writeLines(status, file.path(qa, "status.txt"))
cat(status, "checks=", nrow(results), "\n")
if (!all(results$pass)) quit(save = "no", status = 1L)
