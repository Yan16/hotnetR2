# Analysis-local provenance and logging helpers.

manifest_file_info <- function(path, checksum = FALSE) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  is_file <- exists && !isTRUE(info$isdir)
  md5 <- if (is_file && isTRUE(checksum)) unname(tools::md5sum(path)) else NA_character_
  data.frame(
    path = path,
    exists = exists,
    size_bytes = if (exists) as.numeric(info$size) else NA_real_,
    modified = if (exists) format(info$mtime, "%Y-%m-%d %H:%M:%S %Z") else NA_character_,
    md5 = md5,
    stringsAsFactors = FALSE
  )
}

#' Record an analysis-local workflow event
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param event Short event identifier.
#' @param detail Optional human-readable detail.
#' @param level Event level, normally `"INFO"`, `"WARN"`, or `"ERROR"`.
#' @return Path to `logs/workflow_events.tsv`, invisibly.
#' @export
log_analysis_event <- function(config, event, detail = "", level = "INFO") {
  assert_scalar(event, "event", "character")
  assert_scalar(detail, "detail", "character")
  assert_scalar(level, "level", "character")
  file <- file.path(analysis_paths(config)[["logs"]], "workflow_events.tsv")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    level = level,
    event = event,
    detail = detail,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    row, file = file, sep = "\t", row.names = FALSE, quote = FALSE,
    col.names = !file.exists(file), append = file.exists(file), na = ""
  )
  invisible(file)
}

#' Write a manifest, effective configuration, and R session information
#'
#' All output is restricted to the analysis directory. Checksums are optional:
#' they can be expensive for multi-gigabyte GWAS/reference files, so they are
#' omitted by default but can be requested for a release-quality provenance
#' record.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param checksum Whether to calculate MD5 checksums for existing input files.
#' @return Named paths of written provenance files.
#' @export
write_analysis_provenance <- function(config, checksum = FALSE) {
  assert_flag(checksum, "checksum")
  paths <- analysis_paths(config)
  dir.create(paths[["metadata"]], recursive = TRUE, showWarnings = FALSE)
  inputs <- c(analysis_yaml = config$config_file, input_paths(config))
  manifest <- dplyr::bind_rows(lapply(names(inputs), function(item) {
    cbind(item = item, manifest_file_info(unname(inputs[[item]]), checksum = checksum))
  }))
  manifest_file <- file.path(paths[["metadata"]], "input_manifest.tsv")
  write_tsv_base(manifest, manifest_file)

  session_file <- file.path(paths[["metadata"]], "session_info.txt")
  writeLines(utils::capture.output(utils::sessionInfo()), session_file)
  runtime_file <- file.path(
    paths[["metadata"]],
    paste0("runtime_", config$execution_profile %||% "legacy", ".yaml")
  )
  runtime <- list(
    profile = config$execution_profile %||% "legacy",
    os = unname(Sys.info()[["sysname"]]),
    architecture = R.version$arch,
    r_version = as.character(getRversion()),
    hotnetR2_version = as.character(utils::packageVersion("hotnetR2")),
    selected_execution = config$selected_execution %||% NULL
  )
  ldak_path <- resolve_config_path(config$reference$ldak_exe, config$project_root)
  runtime$ldak <- list(
    configured_path = config$reference$ldak_exe,
    resolved_path = ldak_path,
    symlink_target = Sys.readlink(ldak_path),
    sha256 = if (file.exists(ldak_path) && !isTRUE(file.info(ldak_path)$isdir)) {
      digest::digest(file = ldak_path, algo = "sha256")
    } else {
      NA_character_
    }
  )
  if (file.exists(ldak_path) && file.access(ldak_path, mode = 1L) == 0L) {
    ldak_probe <- probe_executable_version(ldak_path)
    runtime$ldak$observed_version <- ldak_probe$version
    runtime$ldak$version_output <- ldak_probe$output
    runtime$ldak$version_probe_status <- ldak_probe$status
  }
  if (identical(config$hhotnet$execution_mode, "apptainer")) {
    image <- if (is.null(config$hhotnet$apptainer_image)) NA_character_ else
      resolve_config_path(config$hhotnet$apptainer_image, config$project_root)
    runtime$hhotnet <- list(
      mode = "apptainer", resolved_path = image, symlink_target = Sys.readlink(image),
      sha256 = if (!is.na(image) && file.exists(image) && !isTRUE(file.info(image)$isdir)) {
        digest::digest(file = image, algo = "sha256")
      } else {
        NA_character_
      }
    )
  } else {
    source <- if (is.null(config$hhotnet$local_dir)) NA_character_ else
      resolve_config_path(config$hhotnet$local_dir, config$project_root)
    compat <- check_project_compat(
      config$project_root,
      config$requirements$hhotnet_compat %||% ">= 1.0.0"
    )
    runtime$hhotnet <- list(
      mode = "local_python", resolved_path = source,
      symlink_target = Sys.readlink(source), compatibility = compat
    )
  }
  yaml::write_yaml(runtime, runtime_file)
  effective_file <- write_effective_config(config)
  event_file <- log_analysis_event(
    config,
    event = "provenance_written",
    detail = paste0("checksum=", checksum, "; inputs=", nrow(manifest))
  )
  c(
    effective_config = effective_file,
    input_manifest = manifest_file,
    session_info = session_file,
    runtime = runtime_file,
    event_log = event_file
  )
}
