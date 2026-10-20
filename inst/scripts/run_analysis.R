#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: run_analysis.R ANALYSIS_YAML PROFILE STAGE MINIMUM_HOTNETR2_VERSION",
    call. = FALSE
  )
}

config_file <- normalizePath(args[[1L]], mustWork = TRUE)
profile <- args[[2L]]
stage <- args[[3L]]
minimum_version <- args[[4L]]

if (!requireNamespace("hotnetR2", quietly = TRUE)) {
  stop("hotnetR2 is not installed for this R interpreter.", call. = FALSE)
}
installed_version <- as.character(utils::packageVersion("hotnetR2"))
if (utils::compareVersion(installed_version, minimum_version) < 0L) {
  stop(
    "hotnetR2 >= ", minimum_version, " is required, but found ",
    installed_version, ".",
    call. = FALSE
  )
}

config <- hotnetR2::read_analysis_config(config_file, profile = profile)
profile_check <- hotnetR2::check_execution_profile(config)
print(profile_check, row.names = FALSE)
if (!all(profile_check$status == "OK")) {
  stop("Execution-profile preflight failed.", call. = FALSE)
}
inputs <- hotnetR2::check_analysis_inputs(config, write_report = FALSE)
print(inputs, row.names = FALSE)
if (!all(inputs$status == "OK")) {
  stop(
    "Missing required inputs: ",
    paste(inputs$path[inputs$status != "OK"], collapse = ", "),
    call. = FALSE
  )
}

if (identical(config$hhotnet$execution_mode, "local_python")) {
  requirement <- config$requirements$hhotnet_compat
  if (is.null(requirement)) requirement <- ">= 1.0.0"
  compat <- hotnetR2::check_project_compat(config$project_root, requirement)
  print(compat, row.names = FALSE)
  if (!identical(compat$status, "ready")) {
    stop(
      "Shared h-HotNet compatibility installation is ", compat$status,
      ": ", compat$detail,
      call. = FALSE
    )
  }
}

if (identical(stage, "check")) {
  print(hotnetR2::setup_analysis(config, dry_run = TRUE))
  message("PASS: package, schema-v2 YAML, profile, inputs and output boundaries are valid.")
  message("No analysis was run.")
  quit(status = 0L)
}

if (identical(stage, "classification") &&
    !isTRUE(config$enhancer_classification$enabled)) {
  stop("Enhancer classification is disabled in analysis.yaml.", call. = FALSE)
}

paths <- hotnetR2::analysis_paths(config)
dir.create(paths[["metadata"]], recursive = TRUE, showWarnings = FALSE)
writeLines(
  capture.output(sessionInfo()),
  file.path(paths[["metadata"]], paste0("launcher_session_info_", profile, ".txt"))
)

if (identical(stage, "generate-ldak")) {
  hotnetR2::setup_analysis(config, dry_run = FALSE)
  hotnetR2::prepare_ldak_annotations(config, dry_run = FALSE)
  print(hotnetR2::generate_ldak_jobs(config, dry_run = FALSE))
  quit(status = 0L)
}

if (identical(stage, "generate-hotnet")) {
  print(hotnetR2::generate_hhotnet_job(config, dry_run = FALSE))
  quit(status = 0L)
}

if (identical(stage, "all")) {
  stages <- "annotations"
  if (isTRUE(config$enhancer_classification$enabled)) {
    stages <- c(stages, "classification")
  }
  stages <- c(stages, "ldak", "scores", "networks", "hotnet", "exports")
  hotnetR2::run_analysis(config, stages = stages)
} else {
  hotnetR2::run_analysis(config, stages = stage)
}
