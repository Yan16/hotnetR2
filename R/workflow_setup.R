# Setup and input checking for the analysis1-local development package.

write_tsv_base <- function(data, file) {
  utils::write.table(
    data, file = file, sep = "\t", row.names = FALSE, quote = FALSE,
    na = ""
  )
}

#' Check externally provided inputs without modifying them
#'
#' The check reads only the files specified by `analysis.yaml`. Its optional
#' report is written beneath the analysis directory.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param write_report Whether to write `metadata/input_check.tsv`.
#' @return A data frame with one row per required input.
#' @export
check_analysis_inputs <- function(config, write_report = TRUE) {
  inputs <- input_paths(config)
  report <- data.frame(
    item = names(inputs),
    path = unname(inputs),
    exists = file.exists(inputs),
    stringsAsFactors = FALSE
  )
  report$status <- ifelse(report$exists, "OK", "FAIL")
  report$detail <- ifelse(report$exists, "", "missing")

  if (isTRUE(write_report)) {
    paths <- analysis_paths(config)
    dir.create(paths[["metadata"]], recursive = TRUE, showWarnings = FALSE)
    write_tsv_base(report, file.path(paths[["metadata"]], "input_check.tsv"))
    log_analysis_event(
      config,
      event = "input_check",
      detail = paste0(sum(report$status == "OK"), "/", nrow(report), " inputs available"),
      level = if (any(report$status == "FAIL")) "WARN" else "INFO"
    )
  }
  report
}

#' Create the analysis-local workflow directory structure
#'
#' This function can only create directories under the directory holding
#' `analysis.yaml`. It never creates or modifies paths below `project.root`.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param dry_run If `TRUE`, return planned directories without writing them.
#' @param checksum Whether to calculate MD5 checksums when writing provenance.
#' @return Named output paths, invisibly when directories are created.
#' @export
setup_analysis <- function(config, dry_run = TRUE, checksum = FALSE) {
  paths <- analysis_paths(config)
  if (!isTRUE(dry_run)) {
    invisible(lapply(unname(paths[names(paths) != "root"]), dir.create,
                     recursive = TRUE, showWarnings = FALSE))
    write_analysis_provenance(config, checksum = checksum)
    log_analysis_event(config, event = "analysis_setup", detail = "analysis-local directories created")
  }
  paths
}
