# Pre-run validation for h-HotNet input files.

#' Validate selected h-HotNet input files
#'
#' Ensures index IDs are contiguous and unique, every edge endpoint occurs in
#' its index, and the configured score file has unique finite non-negative
#' scores for at least one node in each selected network.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for the configured default, one network name, or a
#'   vector of network names.
#' @return One validation row per selected network; the same table is written
#'   below `hHotnet/summary`.
#' @export
validate_hhotnet_inputs <- function(config, networks = NULL) {
  if (!requireNamespace("readr", quietly = TRUE)) stop("h-HotNet input validation requires readr", call. = FALSE)
  requested <- resolve_networks(config, networks)
  score_name <- config$hhotnet$score_name %||% "nodes_all"
  score_file <- network_data_file(config, paste0(score_name, ".tsv"))
  if (!file.exists(score_file)) stop("Missing h-HotNet score file: ", score_file, call. = FALSE)
  scores <- readr::read_tsv(score_file, col_names = c("gene", "score"),
                            col_types = readr::cols(gene = readr::col_character(), score = readr::col_double()),
                            show_col_types = FALSE, progress = FALSE)
  if (anyDuplicated(scores$gene) || any(!is.finite(scores$score)) || any(scores$score < 0)) {
    stop("h-HotNet score file must have unique finite non-negative scores", call. = FALSE)
  }
  rows <- lapply(requested, function(network) {
    index_file <- network_data_file(config, paste0(network, "_index_gene.tsv"))
    edge_file <- network_data_file(config, paste0(network, "_edge_list.tsv"))
    if (!file.exists(index_file) || !file.exists(edge_file)) stop("Missing h-HotNet input for ", network, call. = FALSE)
    index <- readr::read_tsv(index_file, col_names = c("id", "gene"),
                             col_types = readr::cols(id = readr::col_integer(), gene = readr::col_character()),
                             show_col_types = FALSE, progress = FALSE)
    edges <- readr::read_tsv(edge_file, col_names = c("id1", "id2"),
                             col_types = readr::cols(id1 = readr::col_integer(), id2 = readr::col_integer()),
                             show_col_types = FALSE, progress = FALSE)
    contiguous_ids <- !anyNA(index$id) && identical(index$id, seq_len(nrow(index)))
    unique_genes <- !anyDuplicated(index$gene) && all(!is.na(index$gene) & nzchar(index$gene))
    valid_endpoints <- !anyNA(edges$id1) && !anyNA(edges$id2) && all(edges$id1 %in% index$id) && all(edges$id2 %in% index$id)
    data.frame(
      network = network, index_rows = nrow(index), edge_rows = nrow(edges), score_rows = nrow(scores),
      scored_network_nodes = sum(index$gene %in% scores$gene), contiguous_ids = contiguous_ids,
      unique_genes = unique_genes, valid_endpoints = valid_endpoints,
      ready = contiguous_ids && unique_genes && valid_endpoints && any(index$gene %in% scores$gene),
      stringsAsFactors = FALSE
    )
  })
  report <- do.call(rbind, rows)
  if (any(!report$ready)) stop("One or more selected h-HotNet inputs failed pre-run validation", call. = FALSE)
  output <- network_summary_file(config, "hhotnet_input_validation.tsv")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(report, output)
  report
}
