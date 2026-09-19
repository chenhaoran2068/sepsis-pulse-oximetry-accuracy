s4v_supported_variables <- function() {
  c("spo2", "sao2", "heart_rate", "map", "respiratory_rate", "temperature_c",
    "pao2_mmhg", "paco2_mmhg", "ph", "lactate_mmol_l", "hb_g_dl")
}

s4v_rules <- function() {
  data.frame(
    variable = s4v_supported_variables(),
    lower = c(70, 70, 0, 0, 0, 25, 0, 0, 6.3, 0, 0),
    lower_inclusive = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
    upper = c(100, 100, 300, 300, 100, 45, 1000, 1000, 8.0, Inf, 30),
    upper_inclusive = rep(TRUE, 11L),
    review_lower_exclusive = c(NA, NA, NA, NA, 70, NA, 800, 250, NA, 30, NA),
    review_code = c(NA, NA, NA, NA, "high_respiratory_rate_review",
      "severe_hypothermia_review", "very_high_pao2_review",
      "very_high_paco2_review", NA, "very_high_lactate_review", NA),
    canonical_unit = c("percent", "percent", "beats/min", "mmHg", "breaths/min", "degC",
      "mmHg", "mmHg", "pH scale", "mmol/L", "g/dL"),
    stringsAsFactors = FALSE
  )
}

s4v_validate <- function(value, variable, existing_state = NULL) {
  if (length(variable) != 1L || !variable %in% s4v_supported_variables()) {
    stop("S4V_UNSUPPORTED_VARIABLE", call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(value))
  n <- length(x)
  if (is.null(existing_state)) existing_state <- rep(NA_character_, n)
  if (length(existing_state) != n) stop("S4V_STATE_LENGTH_MISMATCH", call. = FALSE)
  existing_state <- as.character(existing_state)
  state <- existing_state
  state[is.na(x) & (is.na(state) | !nzchar(state))] <- "missing_or_unavailable"
  review <- rep("none", n)

  rule <- s4v_rules()[s4v_rules()$variable == variable, , drop = FALSE]
  nonfinite <- !is.na(x) & !is.finite(x)
  below <- !is.na(x) & is.finite(x) & if (rule$lower_inclusive) x < rule$lower else x <= rule$lower
  above <- !is.na(x) & is.finite(x) & if (rule$upper_inclusive) x > rule$upper else x >= rule$upper
  valid <- !is.na(x) & is.finite(x) & !below & !above

  state[nonfinite] <- "invalid_nonfinite"
  state[below] <- if (rule$lower == 0 && !rule$lower_inclusive) "invalid_nonpositive" else "invalid_below_lower_bound"
  state[above] <- "invalid_above_upper_bound"
  state[valid & (is.na(state) | !nzchar(state))] <- "eligible"
  x[nonfinite | below | above] <- NA_real_

  if (!is.na(rule$review_lower_exclusive) && !is.na(rule$review_code)) {
    review[valid & x > rule$review_lower_exclusive] <- rule$review_code
  }
  if (variable == "temperature_c") {
    review[valid & x <= 30] <- "severe_hypothermia_review"
  }

  data.frame(value = x, value_state = state, diagnostic_review = review,
    stringsAsFactors = FALSE)
}

s4v_apply <- function(data, value_field, variable, state_field = paste0(value_field, "_state"),
                      review_field = paste0(value_field, "_diagnostic_review")) {
  if (!value_field %in% names(data)) stop("S4V_VALUE_FIELD_MISSING", call. = FALSE)
  prior <- if (state_field %in% names(data)) data[[state_field]] else NULL
  checked <- s4v_validate(data[[value_field]], variable, prior)
  data[[value_field]] <- checked$value
  data[[state_field]] <- checked$value_state
  data[[review_field]] <- checked$diagnostic_review
  data
}

s4v_convert_pressure <- function(value, unit) {
  x <- suppressWarnings(as.numeric(value))
  u <- trimws(tolower(as.character(unit)))
  out <- rep(NA_real_, length(x))
  out[u %in% c("mmhg", "mm hg")] <- x[u %in% c("mmhg", "mm hg")]
  out[u == "kpa"] <- x[u == "kpa"] * 7.50062
  out
}

s4v_convert_hb <- function(value, unit) {
  x <- suppressWarnings(as.numeric(value))
  u <- trimws(tolower(as.character(unit)))
  out <- rep(NA_real_, length(x))
  out[u %in% c("g/dl", "g dl-1")] <- x[u %in% c("g/dl", "g dl-1")]
  out[u %in% c("g/l", "g l-1")] <- x[u %in% c("g/l", "g l-1")] / 10
  out
}

s4v_convert_temperature <- function(value, unit) {
  x <- suppressWarnings(as.numeric(value))
  u <- trimws(tolower(as.character(unit)))
  out <- rep(NA_real_, length(x))
  out[u %in% c("c", "degc", "degree c", "degrees celsius")] <- x[u %in% c("c", "degc", "degree c", "degrees celsius")]
  out[u %in% c("f", "degf", "degree f", "degrees fahrenheit")] <- (x[u %in% c("f", "degf", "degree f", "degrees fahrenheit")] - 32) * 5 / 9
  out
}

s4v_assert_no_row_loss <- function(before, after) {
  if (nrow(before) != nrow(after)) stop("S4V_ROW_CARDINALITY_CHANGED", call. = FALSE)
  invisible(TRUE)
}
