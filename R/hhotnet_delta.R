# Manual-delta execution is deliberately separate from permutation completion.
delta_python <- function(config, script, args) {
  if (!identical(config$hhotnet$execution_mode, "local_python")) {
    stop("Manual delta extraction currently requires local_python execution", call. = FALSE)
  }
  compat <- system.file("compat", package = "hotnetR2", mustWork = TRUE)
  temp <- withr::local_tempdir()
  processx::run(config$hhotnet$python_executable %||% "python3", c(script, args),
    env = c("current", PYTHONPATH = paste(c(compat, Sys.getenv("PYTHONPATH")), collapse = .Platform$path.sep),
            MPLCONFIGDIR = temp, XDG_CACHE_HOME = temp), error_on_status = TRUE)
}

run_hhotnet_similarity <- function(config, networks = NULL) {
  networks <- resolve_networks(config, networks)
  ready <- validate_hhotnet_inputs(config, networks)
  if (!all(ready$ready)) stop("HotNet inputs are incomplete", call. = FALSE)
  paths <- analysis_paths(config)
  runner <- generate_hhotnet_job(config, networks, dry_run = FALSE)
  source <- resolve_config_path(config$hhotnet$local_dir, config$project_root)
  tools <- list.files(file.path(source, "src"), "[.]py$", full.names = TRUE)
  reports <- lapply(networks, function(network) {
    folder <- file.path(paths[["hhotnet_intermediate"]], network)
    matrix <- file.path(folder, "similarity_matrix.h5")
    stamp <- file.path(folder, "similarity_fingerprint.rds")
    inputs <- file.path(paths[["hhotnet_data"]], paste0(network, c("_edge_list.tsv", "_index_gene.tsv")))
    signature <- execution_signature(config$hhotnet, c(inputs, tools))
    if (!execution_current(stamp, signature, matrix)) {
      if (dir.exists(folder)) archive_execution_dirs(config, folder)
      processx::run("bash", c(runner, network), error_on_status = TRUE)
      if (!file.exists(matrix) || file.info(matrix)$size == 0) stop("Missing similarity matrix: ", matrix)
      write_execution_stamp(stamp, signature, matrix)
    }
    tibble::tibble(network = network, mode = "similarity_only", similarity_matrix = matrix, complete = TRUE)
  })
  report <- dplyr::bind_rows(reports)
  output <- network_summary_file(config, "hhotnet_similarity_completion.tsv")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(report, output)
  invisible(report)
}

#' Extract exploratory HotNet clusters at manual delta thresholds
#'
#' Requires existing network similarity matrices. Builds only the observed,
#' score-weighted hierarchy once per network (with fingerprinted reuse), then
#' calls the installed HHN `process_hierarchies.cut_hierarchy` function for each
#' delta. No score permutations, null hierarchies, or significance tests are run.
#' Larger delta retains merges at greater heights and usually produces smaller
#' clusters. Heights equal to delta are included, following HHN semantics.
#'
#' @param config Configuration returned by [read_analysis_config()].
#' @param deltas Positive finite thresholds; defaults to 0.05, 0.1, 0.2, 0.5.
#' @param networks Optional configured network names; NULL uses YAML selection.
#' @param min_cluster_size Minimum size included in CX2 and its node/edge tables.
#'   All clusters, including singletons, remain in the raw cluster TSV.
#' @return A manifest with one row per network/delta, paths, counts, and method.
#'   Written under `hHotnet/summary/manual_delta_manifest.tsv`. CX2 contains only
#'   cluster members and within-cluster original network edges, not surrounding
#'   connected components. Outputs are exploratory, not statistically significant.
#' @export
extract_hhotnet_delta_clusters <- function(config, deltas = c(0.05, 0.1, 0.2, 0.5),
                                           networks = NULL, min_cluster_size = 2L) {
  if (!identical(config$hhotnet$execution_mode, "local_python")) {
    stop("Manual delta extraction currently requires local_python execution", call. = FALSE)
  }
  if (!is.numeric(deltas) || !length(deltas) || any(!is.finite(deltas) | deltas <= 0)) {
    stop("deltas must be positive finite numeric values", call. = FALSE)
  }
  if (length(min_cluster_size) != 1L || !is.finite(min_cluster_size) ||
      min_cluster_size < 1 || min_cluster_size != floor(min_cluster_size)) {
    stop("min_cluster_size must be a positive integer", call. = FALSE)
  }
  networks <- resolve_networks(config, networks)
  deltas <- unique(deltas)
  labels <- vapply(deltas, function(x) format(x, scientific = FALSE, trim = TRUE, digits = 15), character(1))
  if (anyDuplicated(labels)) stop("Delta labels collide at 15-digit precision", call. = FALSE)
  paths <- analysis_paths(config)
  source <- file.path(resolve_config_path(config$hhotnet$local_dir, config$project_root), "src")
  helper <- system.file("scripts", "cut_hhotnet_delta.py", package = "hotnetR2", mustWork = TRUE)
  tools <- list.files(source, "[.]py$", full.names = TRUE)
  manifests <- list()
  for (network in networks) {
    matrix <- file.path(paths[["hhotnet_intermediate"]], network, "similarity_matrix.h5")
    if (!file.exists(matrix)) stop("Run the hotnet similarity stage first: ", matrix, call. = FALSE)
    inputs <- file.path(paths[["hhotnet_data"]], c(paste0(network, "_index_gene.tsv"),
                                                paste0(hhotnet_score_name(config), ".tsv")))
    folder <- file.path(paths[["hhotnet_intermediate"]], paste0(network, "_", hhotnet_score_name(config), "_manual_delta"))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    hierarchy <- file.path(folder, "hierarchy_edge_list_0.tsv")
    index <- file.path(folder, "hierarchy_index_gene_0.tsv")
    stamp <- file.path(folder, "hierarchy_fingerprint.rds")
    signature <- execution_signature(config$hhotnet, c(matrix, inputs, tools))
    if (!execution_current(stamp, signature, c(hierarchy, index))) {
      delta_python(config, file.path(source, "construct_hierarchy.py"),
        c("-smf", matrix, "-igf", inputs[[1]], "-gsf", inputs[[2]], "-helf", hierarchy, "-higf", index))
      write_execution_stamp(stamp, signature, c(hierarchy, index))
    }
    result_dir <- file.path(paths[["hhotnet_results"]], network, "manual_delta")
    export_dir <- file.path(paths[["hhotnet_cytoscape"]], "manual_delta")
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(deltas)) {
      suffix <- paste0(network, "_", hhotnet_score_name(config), "_delta", labels[[i]])
      clusters <- file.path(result_dir, paste0("clusters_", suffix, ".tsv"))
      delta_python(config, helper, c("--source", source, "--edges", hierarchy, "--index", index,
                                    "--delta", as.character(deltas[[i]]), "--output", clusters))
      result <- hhotnet_result_data(config, network, cluster_file = clusters)
      nodes <- result$nodes |>
        dplyr::filter(.data$cluster_size >= min_cluster_size) |>
        style_hhotnet_nodes() |>
        dplyr::mutate(delta = deltas[[i]], method = "manual_delta_no_permutations")
      edges <- result$cluster_edges |>
        dplyr::filter(.data$from %in% nodes$name, .data$to %in% nodes$name) |>
        style_hhotnet_edges() |>
        dplyr::mutate(delta = deltas[[i]], method = "manual_delta_no_permutations")
      prefix <- file.path(export_dir, paste0(config$analysis$name, "_", suffix))
      readr::write_tsv(nodes, paste0(prefix, "_nodes.tsv"))
      readr::write_tsv(edges, paste0(prefix, "_edges.tsv"))
      write_hhotnet_cx2(nodes, edges, paste0(prefix, ".cx2"), config$analysis$name, network, TRUE)
      manifests[[length(manifests) + 1L]] <- tibble::tibble(
        network = network, delta = deltas[[i]], method = "manual_delta_no_permutations",
        significance = "not_evaluated", min_cluster_size = min_cluster_size,
        clusters = nrow(result$clusters), exported_clusters = dplyr::n_distinct(nodes$cluster),
        largest_cluster = max(result$clusters$cluster_size), nodes = nrow(nodes), edges = nrow(edges),
        cluster_file = clusters, node_file = paste0(prefix, "_nodes.tsv"),
        edge_file = paste0(prefix, "_edges.tsv"), cx2_file = paste0(prefix, ".cx2"))
    }
  }
  manifest <- dplyr::bind_rows(manifests)
  output <- network_summary_file(config, "manual_delta_manifest.tsv")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(manifest, output)
  manifest
}
