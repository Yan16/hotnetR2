#' Execute regional LDAK with content-aware reuse
#'
#' Reuses completed results only when input, configuration, implementation,
#' and output fingerprints agree. Changed results are archived before a fresh
#' calculation. A symbolic-link result directory is never used as a write target.
#'
#' @param config Configuration returned by [read_analysis_config()].
#' @param node_types Which regional calculations to execute.
#' @return Invisibly, the LDAK structural validation table.
#' @export
run_ldak <- function(config, node_types = c("enhancer", "promoter")) {
  node_types <- match.arg(node_types, c("enhancer", "promoter"), several.ok = TRUE)
  generate_ldak_jobs(config, dry_run = FALSE)
  generate_ldak_local_runner(config, dry_run = FALSE)
  paths <- analysis_paths(config)
  inputs <- input_paths(config)
  inputs <- inputs[!names(inputs) %in% c(
    "aracne", "stringdb_aliases", "stringdb_info", "stringdb_links"
  )]
  for (kind in node_types) {
    result <- file.path(paths[["ldak_results"]], ldak_result_subdir(config, kind))
    assert_unlinked_result(result)
    signature <- execution_signature(config$ldak, c(inputs,
      file.path(paths[["ldak_annotations"]], paste0(if (kind == "enhancer") "enhancers" else "promoters", ".loc"))))
    stamp <- file.path(paths[["metadata"]], paste0("ldak_", kind, "_fingerprint.rds"))
    artifacts <- file.path(result, c("remls.all", "genes.details"))
    if (execution_current(stamp, signature, artifacts)) next
    archive_execution_dirs(config, result)
    log <- file.path(paths[["logs"]], paste0("ldak_", kind, ".log"))
    processx::run("bash", file.path(paths[["ldak_scripts"]], paste0("run_ldak_", kind, ".sh")),
                  stdout = log, stderr = "2>&1", error_on_status = TRUE)
    # LDAK can report input errors without a nonzero exit code.
    if (!all(file.exists(artifacts)) || any(file.info(artifacts)$size == 0)) {
      stop("LDAK did not produce complete results; inspect ", log, call. = FALSE)
    }
    write_execution_stamp(stamp, signature, artifacts)
  }
  invisible(validate_ldak_results(config, node_types))
}

assert_unlinked_result <- function(path) {
  target <- Sys.readlink(path)
  if (!is.na(target) && nzchar(target)) {
    stop("Refusing to write through linked LDAK results: ", path, call. = FALSE)
  }
}

#' Execute configured hierarchical HotNet networks with fingerprinted resume
#'
#' @param config Configuration returned by [read_analysis_config()].
#' @param networks Optional configured network names; NULL selects YAML defaults.
#' @return Invisibly, the completion validation table.
#' @export
run_hhotnet <- function(config, networks = NULL) {
  if (identical(config$hhotnet$analysis_mode, "manual_delta")) {
    return(run_hhotnet_similarity(config, networks))
  }
  networks <- resolve_networks(config, networks)
  ready <- validate_hhotnet_inputs(config, networks)
  if (!all(ready$ready)) stop("HotNet inputs are incomplete", call. = FALSE)
  runner <- generate_hhotnet_job(config, networks, dry_run = FALSE)
  paths <- analysis_paths(config)
  if (identical(config$hhotnet$execution_mode, "local_python")) {
    source <- resolve_config_path(config$hhotnet$local_dir, config$project_root)
    compat <- resolve_config_path(config$hhotnet$compat_dir %||%
      file.path(paths[["hhotnet"]], "compat"), config$project_root)
    tool_files <- c(list.files(file.path(source, "src"), "[.]py$", full.names = TRUE, recursive = TRUE),
                    list.files(compat, "[.]py$", full.names = TRUE, recursive = TRUE))
    environment <- processx::run(config$hhotnet$python_executable %||% "python3",
      c("-c", "import sys, importlib.metadata as m; print(sys.version); print(sorted((d.metadata['Name'], d.version) for d in m.distributions()))"))$stdout
  } else {
    tool_files <- resolve_config_path(config$hhotnet$apptainer_image, config$project_root)
    environment <- processx::run("apptainer", "--version")$stdout
  }
  # Generated shell runners support overrides; include them in reuse identity.
  overrides <- Sys.getenv()
  overrides <- overrides[grepl("^(HHOTNET_|SLURM_CPUS_PER_TASK$|PYTHONPATH$)", names(overrides))]
  for (network in networks) {
    inputs <- file.path(paths[["hhotnet_data"]],
                         c(paste0(network, "_index_gene.tsv"), paste0(network, "_edge_list.tsv"),
                           paste0(config$hhotnet$score_name, ".tsv")))
    signature <- execution_signature(list(config$hhotnet, environment, overrides), c(inputs, tool_files))
    stamp <- file.path(paths[["metadata"]], paste0(network, "_fingerprint.rds"))
    result <- file.path(paths[["hhotnet_results"]], network)
    files <- list.files(result, full.names = TRUE, recursive = TRUE)
    if (length(files) && execution_current(stamp, signature, files)) next
    archive_execution_dirs(config, c(result,
      file.path(paths[["hhotnet_intermediate"]], network),
      file.path(paths[["hhotnet_intermediate"]], paste0(network, "_", config$hhotnet$score_name))))
    processx::run("bash", c(runner, network), error_on_status = TRUE)
    complete <- validate_hhotnet_completion(config, network)
    if (!all(complete$complete)) stop("HotNet did not complete ", network, call. = FALSE)
    write_execution_stamp(stamp, signature, list.files(result, full.names = TRUE, recursive = TRUE))
  }
  invisible(validate_hhotnet_completion(config, networks))
}

execution_signature <- function(settings, files) {
  implementation <- ls(envir = asNamespace("hotnetR2"), all.names = TRUE)
  implementation <- purrr::map(implementation, ~ get(.x, envir = asNamespace("hotnetR2")))
  implementation <- purrr::keep(implementation, is.function)
  digest::digest(list(settings, purrr::map_chr(files, ~ digest::digest(file = .x, algo = "sha256")),
                      purrr::map(implementation, body)), algo = "sha256")
}

execution_current <- function(stamp, signature, files) {
  if (!file.exists(stamp) || !all(file.exists(files))) return(FALSE)
  previous <- readRDS(stamp)
  identical(previous$signature, signature) &&
    identical(previous$outputs, purrr::map_chr(files, ~ digest::digest(file = .x, algo = "sha256")))
}

write_execution_stamp <- function(stamp, signature, files) {
  dir.create(dirname(stamp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(signature = signature,
                outputs = purrr::map_chr(files, ~ digest::digest(file = .x, algo = "sha256"))), stamp)
}

archive_execution_dirs <- function(config, dirs) {
  dirs <- dirs[dir.exists(dirs)]
  if (!length(dirs)) return(invisible(NULL))
  archive <- tempfile("superseded-", tmpdir = analysis_paths(config)[["root"]])
  dir.create(archive)
  for (i in seq_along(dirs)) {
    destination <- file.path(archive, paste0(i, "-", basename(dirs[[i]])))
    if (!file.rename(dirs[[i]], destination)) stop("Cannot archive ", dirs[[i]], call. = FALSE)
  }
  message("Archived superseded results in ", archive)
}

#' Run the installed v2/v3 workflow
#'
#' @param config A configuration list or YAML path.
#' @param stages Ordered stages to execute; defaults to the complete workflow.
#' @return Invisibly, the resolved configuration. Writes analysis outputs.
#' @export
run_analysis <- function(config, stages = c("annotations", "classification", "ldak",
                                            "scores", "networks", "hotnet", "exports")) {
  if (is.character(config)) config <- read_analysis_config(config)
  stages <- match.arg(stages, c("annotations", "classification", "ldak", "scores",
                               "networks", "hotnet", "exports"), several.ok = TRUE)
  setup_analysis(config, dry_run = FALSE)
  for (stage in stages) {
    message("hotnetR2: ", stage)
    switch(stage,
      annotations = prepare_ldak_annotations(config, dry_run = FALSE),
      classification = classify_enhancers(config),
      ldak = run_ldak(config),
      scores = summarize_analysis_ldak_results(config, dry_run = FALSE),
      networks = build_hhotnet_networks(config, dry_run = FALSE),
      hotnet = run_hhotnet(config),
      exports = {
        if (identical(config$hhotnet$analysis_mode, "manual_delta")) {
          extract_hhotnet_delta_clusters(config,
            deltas = unlist(config$hhotnet$delta_values %||% c(0.05, 0.1, 0.2, 0.5)))
          next
        }
        summarize_hhotnet_results(config)
        export_hhotnet_graphs(config)
        valid <- validate_hhotnet_exports(config)
        if (!all(valid$valid)) stop("Graph export validation failed", call. = FALSE)
      }
    )
  }
  invisible(config)
}
