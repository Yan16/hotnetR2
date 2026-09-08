# Completion checks for locally run h-HotNet analyses.

#' Validate selected h-HotNet execution checkpoints and outputs
#'
#' Reports whether every expected permutation hierarchy, cluster table, and
#' size plot exists for each selected network. Incomplete networks remain
#' restartable by the generated runner.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for the configured default, one network name, or a
#'   vector of network names.
#' @return A completion data frame written below `hHotnet/summary`.
#' @export
validate_hhotnet_completion <- function(config, networks = NULL) {
  requested <- resolve_networks(config, networks)
  paths <- analysis_paths(config)
  score_name <- config$hhotnet$score_name %||% "nodes_all"
  expected <- as.integer(config$hhotnet$num_permutations) + 1L
  rows <- lapply(requested, function(network) {
    score_dir <- file.path(paths[["hhotnet_intermediate"]], paste0(network, "_", score_name))
    results_dir <- file.path(paths[["hhotnet_results"]], network)
    hierarchy_edges <- length(list.files(score_dir, pattern = "^hierarchy_edge_list_[0-9]+\\.tsv$"))
    hierarchy_indices <- length(list.files(score_dir, pattern = "^hierarchy_index_gene_[0-9]+\\.tsv$"))
    cluster_file <- file.path(results_dir, paste0("clusters_", network, "_", score_name, ".tsv"))
    plot_file <- file.path(results_dir, paste0("sizes_", network, "_", score_name, ".pdf"))
    data.frame(
      network = network, expected_hierarchies = expected,
      hierarchy_edge_files = hierarchy_edges, hierarchy_index_files = hierarchy_indices,
      cluster_file = cluster_file, cluster_exists = file.exists(cluster_file) && file.info(cluster_file)$size > 0,
      plot_file = plot_file, plot_exists = file.exists(plot_file) && file.info(plot_file)$size > 0,
      complete = hierarchy_edges == expected && hierarchy_indices == expected &&
        file.exists(cluster_file) && file.info(cluster_file)$size > 0 && file.exists(plot_file) && file.info(plot_file)$size > 0,
      stringsAsFactors = FALSE
    )
  })
  report <- do.call(rbind, rows)
  output <- network_summary_file(config, "hhotnet_completion_validation.tsv")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(report, output)
  report
}
