r05_stop <- function(code) stop(code, call. = FALSE)

r05_event_columns <- c(
  "stay_key_internal", "patient_key_internal", "event_time_min", "saturation_percent",
  "source_family", "source_priority", "range_status", "analysis_eligible_range",
  "artifact_qc_status", "possible_transient_artifact", "stable_event_ordinal"
)

r05_empty_pairs <- function() data.frame(
  cohort = character(), patient_key_internal = character(), stay_key_internal = character(),
  spo2_stable_event_ordinal = integer(), sao2_stable_event_ordinal = integer(),
  spo2_relative_icu_minutes = numeric(), sao2_relative_icu_minutes = numeric(),
  spo2_saturation_percent = numeric(), sao2_saturation_percent = numeric(),
  absolute_lag_minutes = numeric(), pair_direction = character(), pair_window_minutes = numeric(),
  contains_possible_transient_artifact = logical(), stringsAsFactors = FALSE
)

r05_empty_component_evidence <- function() data.frame(
  cohort = character(), window_minutes = numeric(), component_ordinal = integer(),
  candidate_edges = double(), spo2_nodes = double(), sao2_nodes = double(),
  pre_artifact_pair_rows = double(), feasibility_calls = double(), elapsed_seconds = double(),
  stringsAsFactors = FALSE
)

r05_validate_events <- function(events, cohort, modality) {
  if (!is.data.frame(events) || length(setdiff(r05_event_columns, names(events)))) r05_stop("R05_EVENT_SCHEMA_ERROR")
  forbidden <- c("stay_raw", "patient_raw", "subject_id", "hadm_id", "stay_id", "comment", "specimen_label", "event_time_absolute")
  if (length(intersect(forbidden, names(events)))) r05_stop("R05_RAW_FIELD_REJECTED")
  text_columns <- c("stay_key_internal", "patient_key_internal", "source_family", "range_status", "artifact_qc_status")
  if (any(vapply(events[text_columns], function(x) !is.character(x), logical(1)))) r05_stop("R05_EVENT_TYPE_ERROR")
  if (!is.numeric(events$event_time_min) || !is.numeric(events$saturation_percent) || !is.numeric(events$stable_event_ordinal)) r05_stop("R05_EVENT_TYPE_ERROR")
  if (!is.logical(events$analysis_eligible_range) || !is.logical(events$possible_transient_artifact)) r05_stop("R05_EVENT_TYPE_ERROR")
  if (any(!nzchar(events$stay_key_internal)) || any(!nzchar(events$patient_key_internal)) ||
      any(!is.finite(events$event_time_min)) || any(!is.finite(events$saturation_percent)) ||
      any(!is.finite(events$stable_event_ordinal)) || any(events$stable_event_ordinal < 1) ||
      any(abs(events$stable_event_ordinal - round(events$stable_event_ordinal)) > 1e-8)) r05_stop("R05_EVENT_VALUE_ERROR")
  if (any(!(events$range_status %in% c("below_70", "eligible_70_100", "above_100")))) r05_stop("R05_RANGE_STATUS_ERROR")
  if (!identical(as.logical(events$analysis_eligible_range), events$range_status == "eligible_70_100")) r05_stop("R05_RANGE_INTERFACE_ERROR")
  if (anyDuplicated(paste(events$stay_key_internal, events$stable_event_ordinal, sep = "\r"))) r05_stop("R05_EVENT_ORDINAL_DUPLICATE")
  mapping <- unique(events[c("stay_key_internal", "patient_key_internal")])
  if (anyDuplicated(mapping$stay_key_internal)) r05_stop("R05_STAY_PATIENT_INTERFACE_ERROR")
  invisible(TRUE)
}

r05_validate_modalities <- function(spo2, sao2) {
  spo_map <- unique(spo2[c("stay_key_internal", "patient_key_internal")])
  sao_map <- unique(sao2[c("stay_key_internal", "patient_key_internal")])
  both <- merge(spo_map, sao_map, by = "stay_key_internal", suffixes = c("_spo", "_sao"))
  if (nrow(both) && any(both$patient_key_internal_spo != both$patient_key_internal_sao)) r05_stop("R05_CROSS_MODALITY_PATIENT_MISMATCH")
  invisible(TRUE)
}

r05_pair_signature <- function(pairs) {
  if (!nrow(pairs)) return("")
  z <- pairs[order(pairs$sao2_stable_event_ordinal, pairs$spo2_stable_event_ordinal), , drop = FALSE]
  paste(paste(sprintf("%020d", z$sao2_stable_event_ordinal), sprintf("%020d", z$spo2_stable_event_ordinal), sep = "~"), collapse = "|")
}

r05_build_edges <- function(spo2, sao2, window_minutes) {
  spo2 <- spo2[order(spo2$event_time_min, spo2$stable_event_ordinal), , drop = FALSE]
  sao2 <- sao2[order(sao2$event_time_min, sao2$stable_event_ordinal), , drop = FALSE]
  if (!nrow(spo2) || !nrow(sao2)) return(list(spo2 = spo2, sao2 = sao2, edges = data.frame(spo_index = integer(), sao_index = integer(), lag_minutes = numeric())))
  blocks <- vector("list", nrow(sao2))
  for (j in seq_len(nrow(sao2))) {
    lower <- findInterval(sao2$event_time_min[[j]] - window_minutes - 1e-10, spo2$event_time_min) + 1L
    upper <- findInterval(sao2$event_time_min[[j]] + window_minutes + 1e-10, spo2$event_time_min)
    if (lower <= upper) {
      indices <- seq.int(lower, upper)
      blocks[[j]] <- data.frame(spo_index = indices, sao_index = rep.int(j, length(indices)),
        lag_minutes = abs(spo2$event_time_min[indices] - sao2$event_time_min[[j]]), stringsAsFactors = FALSE)
    }
  }
  blocks <- Filter(Negate(is.null), blocks)
  edges <- if (length(blocks)) do.call(rbind, blocks) else data.frame(spo_index = integer(), sao_index = integer(), lag_minutes = numeric())
  if (nrow(edges) && any(!is.finite(edges$lag_minutes) | edges$lag_minutes > window_minutes + 1e-8)) r05_stop("R05_EDGE_WINDOW_ERROR")
  list(spo2 = spo2, sao2 = sao2, edges = edges)
}

r05_component_edges <- function(edges) {
  if (!nrow(edges)) return(list())
  graph <- igraph::graph_from_data_frame(data.frame(from = paste0("P", edges$spo_index), to = paste0("A", edges$sao_index), stringsAsFactors = FALSE), directed = FALSE)
  membership <- igraph::components(graph)$membership
  edge_component <- as.integer(membership[paste0("P", edges$spo_index)])
  ids <- sort(unique(edge_component))
  order_ids <- ids[order(vapply(ids, function(id) min(edges$sao_index[edge_component == id]), numeric(1)))]
  lapply(order_ids, function(id) edges[edge_component == id, , drop = FALSE])
}

r05_optimize_edges <- function(edges, n_spo2, n_sao2) {
  if (!nrow(edges)) return(list(pair_count = 0L, total_lag_minutes = 0, selected_edge_ids = integer()))
  max_pairs <- min(n_spo2, n_sao2)
  max_lag <- max(edges$lag_minutes)
  pair_bonus <- as.double(max_pairs + 1L) * (max_lag + 1)
  if (!is.finite(pair_bonus) || pair_bonus <= max_lag) r05_stop("R05_WEIGHT_RANGE_ERROR")
  vertices <- data.frame(name = c(paste0("P", seq_len(n_spo2)), paste0("A", seq_len(n_sao2))),
    type = c(rep(TRUE, n_spo2), rep(FALSE, n_sao2)), stringsAsFactors = FALSE)
  graph <- igraph::graph_from_data_frame(data.frame(from = paste0("P", edges$spo_index), to = paste0("A", edges$sao_index), stringsAsFactors = FALSE), directed = FALSE, vertices = vertices)
  graph <- igraph::set_edge_attr(graph, "weight", value = pair_bonus - edges$lag_minutes)
  solved <- igraph::max_bipartite_match(graph, types = igraph::vertex_attr(graph, "type"), weights = igraph::edge_attr(graph, "weight"))
  matching <- solved$matching[paste0("P", seq_len(n_spo2))]
  matched_spo <- which(!is.na(matching) & nzchar(matching))
  if (!length(matched_spo)) return(list(pair_count = 0L, total_lag_minutes = 0, selected_edge_ids = integer()))
  matched_sao <- as.integer(sub("^A", "", matching[matched_spo]))
  selected <- match(paste(matched_spo, matched_sao, sep = "\r"), paste(edges$spo_index, edges$sao_index, sep = "\r"))
  if (anyNA(selected)) r05_stop("R05_MATCHING_INTERFACE_ERROR")
  list(pair_count = as.integer(length(selected)), total_lag_minutes = sum(edges$lag_minutes[selected]), selected_edge_ids = as.integer(selected))
}

r05_feasible_after_constraints <- function(edges, n_spo2, n_sao2, target_count, target_lag, fixed_edge_ids, excluded_sao_indices) {
  if (length(fixed_edge_ids) && (anyDuplicated(edges$spo_index[fixed_edge_ids]) || anyDuplicated(edges$sao_index[fixed_edge_ids]))) return(FALSE)
  if (length(fixed_edge_ids) && any(edges$sao_index[fixed_edge_ids] %in% excluded_sao_indices)) return(FALSE)
  fixed_count <- length(fixed_edge_ids)
  fixed_lag <- if (fixed_count) sum(edges$lag_minutes[fixed_edge_ids]) else 0
  if (fixed_count > target_count || fixed_lag > target_lag + 1e-8) return(FALSE)
  allowed <- rep(TRUE, nrow(edges))
  if (length(excluded_sao_indices)) allowed <- allowed & !(edges$sao_index %in% excluded_sao_indices)
  if (fixed_count) allowed <- allowed & !(edges$spo_index %in% edges$spo_index[fixed_edge_ids]) & !(edges$sao_index %in% edges$sao_index[fixed_edge_ids])
  solved <- r05_optimize_edges(edges[allowed, , drop = FALSE], n_spo2, n_sao2)
  fixed_count + solved$pair_count == target_count && abs(fixed_lag + solved$total_lag_minutes - target_lag) <= 1e-8
}

r05_emit_pairs <- function(spo2, sao2, edges, selected_edge_ids, cohort, window_minutes) {
  if (!length(selected_edge_ids)) return(r05_empty_pairs())
  ordered <- selected_edge_ids[order(edges$sao_index[selected_edge_ids])]
  si <- edges$spo_index[ordered]; ai <- edges$sao_index[ordered]
  direction <- ifelse(spo2$event_time_min[si] < sao2$event_time_min[ai], "spo2_before_sao2",
    ifelse(spo2$event_time_min[si] > sao2$event_time_min[ai], "spo2_after_sao2", "same_time"))
  data.frame(cohort = cohort, patient_key_internal = spo2$patient_key_internal[si], stay_key_internal = spo2$stay_key_internal[si],
    spo2_stable_event_ordinal = as.integer(spo2$stable_event_ordinal[si]), sao2_stable_event_ordinal = as.integer(sao2$stable_event_ordinal[ai]),
    spo2_relative_icu_minutes = spo2$event_time_min[si], sao2_relative_icu_minutes = sao2$event_time_min[ai],
    spo2_saturation_percent = spo2$saturation_percent[si], sao2_saturation_percent = sao2$saturation_percent[ai],
    absolute_lag_minutes = edges$lag_minutes[ordered], pair_direction = direction, pair_window_minutes = window_minutes,
    contains_possible_transient_artifact = spo2$possible_transient_artifact[si], stringsAsFactors = FALSE)
}

r05_match_component <- function(spo2, sao2, edges, cohort, window_minutes, component_edge_limit, feasibility_call_limit) {
  if (nrow(edges) > component_edge_limit) r05_stop("R05_COMPONENT_EDGE_LIMIT_REJECTED")
  baseline <- r05_optimize_edges(edges, nrow(spo2), nrow(sao2))
  if (!baseline$pair_count) return(list(pairs = r05_empty_pairs(), feasibility_calls = 0L))
  fixed <- integer(); excluded_sao <- integer(); last_sao <- 0L; calls <- 0L
  for (position in seq_len(baseline$pair_count)) {
    candidates <- which(!(edges$sao_index %in% excluded_sao) & edges$sao_index > last_sao &
      !(edges$spo_index %in% edges$spo_index[fixed]) & !(edges$sao_index %in% edges$sao_index[fixed]))
    candidates <- candidates[order(sao2$stable_event_ordinal[edges$sao_index[candidates]], spo2$stable_event_ordinal[edges$spo_index[candidates]])]
    chosen <- NA_integer_
    for (edge_id in candidates) {
      calls <- calls + 1L
      if (calls > feasibility_call_limit) r05_stop("R05_FEASIBILITY_LIMIT_REJECTED")
      prior_sao <- if (edges$sao_index[[edge_id]] > last_sao + 1L) seq.int(last_sao + 1L, edges$sao_index[[edge_id]] - 1L) else integer()
      proposed_excluded <- sort(unique(c(excluded_sao, prior_sao)))
      proposed_fixed <- c(fixed, edge_id)
      if (r05_feasible_after_constraints(edges, nrow(spo2), nrow(sao2), baseline$pair_count, baseline$total_lag_minutes, proposed_fixed, proposed_excluded)) {
        chosen <- edge_id; fixed <- proposed_fixed; excluded_sao <- proposed_excluded; last_sao <- edges$sao_index[[edge_id]]; break
      }
    }
    if (is.na(chosen)) r05_stop("R05_TERTIARY_TIEBREAK_ERROR")
  }
  pairs <- r05_emit_pairs(spo2, sao2, edges, fixed, cohort, window_minutes)
  if (nrow(pairs) != baseline$pair_count || abs(sum(pairs$absolute_lag_minutes) - baseline$total_lag_minutes) > 1e-8) r05_stop("R05_POSTCONDITION_ERROR")
  list(pairs = pairs, feasibility_calls = as.integer(calls))
}

r05_pair_one_stay <- function(spo2, sao2, cohort, window_minutes, component_edge_limit, feasibility_call_limit, decompose = TRUE) {
  built <- r05_build_edges(spo2, sao2, window_minutes)
  if (!nrow(built$edges)) return(list(pairs = r05_empty_pairs(), candidate_edges = 0L, components = r05_empty_component_evidence(), feasibility_calls = 0L))
  edge_sets <- if (decompose) r05_component_edges(built$edges) else list(built$edges)
  selected <- list(); evidence <- list(); calls <- 0L
  for (rank in seq_along(edge_sets)) {
    edge_set <- edge_sets[[rank]]; spi <- sort(unique(edge_set$spo_index)); sai <- sort(unique(edge_set$sao_index))
    component_spo <- built$spo2[spi, , drop = FALSE]; component_sao <- built$sao2[sai, , drop = FALSE]
    remapped <- data.frame(spo_index = match(edge_set$spo_index, spi), sao_index = match(edge_set$sao_index, sai), lag_minutes = edge_set$lag_minutes, stringsAsFactors = FALSE)
    started <- proc.time()[["elapsed"]]
    solved <- r05_match_component(component_spo, component_sao, remapped, cohort, window_minutes, component_edge_limit, feasibility_call_limit)
    elapsed <- proc.time()[["elapsed"]] - started; calls <- calls + solved$feasibility_calls
    evidence[[rank]] <- data.frame(cohort = cohort, window_minutes = window_minutes, component_ordinal = as.integer(rank), candidate_edges = as.double(nrow(remapped)), spo2_nodes = as.double(nrow(component_spo)), sao2_nodes = as.double(nrow(component_sao)), pre_artifact_pair_rows = as.double(nrow(solved$pairs)), feasibility_calls = as.double(solved$feasibility_calls), elapsed_seconds = as.double(elapsed), stringsAsFactors = FALSE)
    if (nrow(solved$pairs)) selected[[length(selected) + 1L]] <- solved$pairs
  }
  list(pairs = if (length(selected)) do.call(rbind, selected) else r05_empty_pairs(), candidate_edges = as.integer(nrow(built$edges)), components = do.call(rbind, evidence), feasibility_calls = as.integer(calls))
}

r05_ranges <- function(events) {
  x <- events[order(events$stay_key_internal, events$event_time_min, events$stable_event_ordinal), , drop = FALSE]
  if (!nrow(x)) return(list(events = x, index = data.frame(stay_key_internal = character(), start_row = integer(), end_row = integer(), stringsAsFactors = FALSE)))
  runs <- rle(x$stay_key_internal)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  list(events = x, index = data.frame(stay_key_internal = runs$values, start_row = starts, end_row = ends, stringsAsFactors = FALSE))
}

r05_pair_cohort <- function(spo2, sao2, cohort, window_minutes, component_edge_limit = 50000L, feasibility_call_limit = 500000L, decompose = TRUE) {
  r05_validate_events(spo2, cohort, "SpO2"); r05_validate_events(sao2, cohort, "SaO2"); r05_validate_modalities(spo2, sao2)
  spo2 <- spo2[spo2$analysis_eligible_range, , drop = FALSE]; sao2 <- sao2[sao2$analysis_eligible_range, , drop = FALSE]
  sp <- r05_ranges(spo2); sa <- r05_ranges(sao2)
  common <- merge(sp$index, sa$index, by = "stay_key_internal", suffixes = c("_spo", "_sao"))
  selected <- list(); evidence <- list(); total_edges <- 0; total_calls <- 0L; component_offset <- 0L
  if (nrow(common)) common <- common[order(common$stay_key_internal), , drop = FALSE]
  for (row in seq_len(nrow(common))) {
    solved <- r05_pair_one_stay(sp$events[common$start_row_spo[[row]]:common$end_row_spo[[row]], , drop = FALSE], sa$events[common$start_row_sao[[row]]:common$end_row_sao[[row]], , drop = FALSE], cohort, window_minutes, component_edge_limit, feasibility_call_limit, decompose)
    total_edges <- total_edges + solved$candidate_edges; total_calls <- total_calls + solved$feasibility_calls
    if (nrow(solved$components)) { solved$components$component_ordinal <- solved$components$component_ordinal + component_offset; component_offset <- max(solved$components$component_ordinal); evidence[[length(evidence) + 1L]] <- solved$components }
    if (nrow(solved$pairs)) selected[[length(selected) + 1L]] <- solved$pairs
  }
  pre <- if (length(selected)) do.call(rbind, selected) else r05_empty_pairs()
  if (nrow(pre) && (anyDuplicated(paste(pre$stay_key_internal, pre$spo2_stable_event_ordinal, sep = "\r")) || anyDuplicated(paste(pre$stay_key_internal, pre$sao2_stable_event_ordinal, sep = "\r")) || any(pre$absolute_lag_minutes > window_minutes + 1e-8))) r05_stop("R05_ONE_TO_ONE_POSTCONDITION_ERROR")
  final <- pre[!pre$contains_possible_transient_artifact, , drop = FALSE]
  aggregate <- data.frame(cohort = cohort, window_minutes = window_minutes, eligible_spo2_events = as.double(nrow(spo2)), eligible_sao2_events = as.double(nrow(sao2)), common_stays = as.double(nrow(common)), candidate_edges = as.double(total_edges), temporal_components = as.double(component_offset), pre_artifact_pair_count = as.double(nrow(pre)), pairs_with_possible_transient_artifact = as.double(sum(pre$contains_possible_transient_artifact)), final_retained_pair_count = as.double(nrow(final)), feasibility_calls = as.double(total_calls), artifact_qc_mode = if (identical(cohort, "SICDB")) "NOT_ASSESSABLE_TEMPORAL_RESOLUTION" else "P99_15_MIN_THREE_POINT", stringsAsFactors = FALSE)
  direction <- data.frame(cohort = cohort, window_minutes = window_minutes, pair_stage = rep(c("pre_artifact", "final"), each = 3L), pair_direction = rep(c("spo2_before_sao2", "same_time", "spo2_after_sao2"), 2L), n = c(vapply(c("spo2_before_sao2", "same_time", "spo2_after_sao2"), function(x) sum(pre$pair_direction == x), numeric(1)), vapply(c("spo2_before_sao2", "same_time", "spo2_after_sao2"), function(x) sum(final$pair_direction == x), numeric(1))), stringsAsFactors = FALSE)
  list(pre_artifact_pairs = pre, final_pairs = final, aggregate = aggregate, direction = direction, component_evidence = if (length(evidence)) do.call(rbind, evidence) else r05_empty_component_evidence())
}
