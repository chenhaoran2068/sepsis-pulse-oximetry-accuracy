args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: check_r_packages.R RUN_LIBRARY RECEIPT_PATH", call. = FALSE)
run_library <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
receipt <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(run_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(run_library, .libPaths()))
needed <- c("data.table", "igraph", "digest", "lubridate", "stringr", "mice",
            "miceadds", "lme4", "Matrix", "nlme", "car", "sandwich",
            "ranger", "fastshap")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  repository <- Sys.getenv("DEMO_CRAN_REPOSITORY", "https://cloud.r-project.org")
  cat("Installing missing R packages into isolated run library:", paste(missing, collapse = ", "), "\n")
  tryCatch(
    utils::install.packages(missing, lib = run_library, repos = repository, quiet = TRUE),
    error = function(e) stop("DEMO_R_DEPENDENCY_INSTALL_FAILED: ", conditionMessage(e), call. = FALSE)
  )
}
still_missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("DEMO_R_DEPENDENCY_MISSING: ", paste(still_missing, collapse = ", "),
       ". Check network access to the configured CRAN repository or install into an R library.", call. = FALSE)
}
versions <- vapply(needed, function(package) as.character(utils::packageVersion(package)), character(1))
write.table(data.frame(package = needed, version = versions), receipt,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("DEMO_R_DEPENDENCIES_PASS\n")
