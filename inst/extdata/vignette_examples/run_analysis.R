args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_analysis.R ANALYSIS_YAML STAGE EXPECTED_VERSION", call. = FALSE)
}

config_file <- normalizePath(args[[1L]], mustWork = TRUE)
stage <- args[[2L]]
expected_version <- args[[3L]]

if (!requireNamespace("hotnetR2", quietly = TRUE)) {
  stop("hotnetR2 is not installed for this R interpreter.", call. = FALSE)
}
installed_version <- as.character(utils::packageVersion("hotnetR2"))
if (!identical(installed_version, expected_version)) {
  stop(
    "Expected hotnetR2 ", expected_version, ", but found ", installed_version,
    ". Set HOTNETR2_VERSION only for an intentional version change.",
    call. = FALSE
  )
}

config <- hotnetR2::read_analysis_config(config_file)
inputs <- hotnetR2::check_analysis_inputs(config, write_report = FALSE)
print(inputs, row.names = FALSE)
if (!all(inputs$status == "OK")) {
  stop(
    "Missing required inputs: ",
    paste(inputs$path[inputs$status != "OK"], collapse = ", "),
    call. = FALSE
  )
}

if (identical(stage, "check")) {
  print(hotnetR2::setup_analysis(config, dry_run = TRUE))
  message("PASS: package, YAML, inputs and output boundaries are valid.")
  message("No analysis was run.")
  quit(status = 0L)
}

if (identical(stage, "classification") &&
    !isTRUE(config$enhancer_classification$enabled)) {
  stop("Enhancer classification is disabled in analysis.yaml.", call. = FALSE)
}

if (identical(config$hhotnet$execution_mode, "local_python")) {
  compat_source <- system.file("compat", package = "hotnetR2", mustWork = TRUE)
  compat_target <- config$hhotnet$compat_dir
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", compat_target)) {
    compat_target <- file.path(config$project_root, compat_target)
  }
  dir.create(compat_target, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(
    list.files(compat_source, full.names = TRUE, all.files = TRUE, no.. = TRUE),
    compat_target, recursive = TRUE, overwrite = TRUE,
    copy.mode = TRUE, copy.date = TRUE
  )
  if (!all(copied)) stop("Could not copy HHN compatibility files.", call. = FALSE)
}

paths <- hotnetR2::analysis_paths(config)
dir.create(paths[["metadata"]], recursive = TRUE, showWarnings = FALSE)
writeLines(
  capture.output(sessionInfo()),
  file.path(paths[["metadata"]], "launcher_session_info.txt")
)

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
