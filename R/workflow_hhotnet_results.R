# Result parsing, summaries, and portable graph exports for completed h-HotNet runs.

require_hhotnet_result_packages <- function(include_json = FALSE) {
  packages <- c("dplyr", "readr", "igraph")
  if (isTRUE(include_json)) packages <- c(packages, "jsonlite")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("h-HotNet result processing requires package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

hhotnet_score_name <- function(config) {
  as.character(config$hhotnet$score_name %||% "nodes_all")
}

hhotnet_cluster_result_file <- function(config, network) {
  file.path(
    analysis_paths(config)[["hhotnet_results"]], network,
    paste0("clusters_", network, "_", hhotnet_score_name(config), ".tsv")
  )
}

parse_hhotnet_cluster_headers <- function(lines, network) {
  lines <- sub("^#\\s*", "", lines)
  matches <- regexec("^([^:]+):\\s*(.*)$", lines)
  parts <- regmatches(lines, matches)
  valid <- lengths(parts) == 3L
  if (!any(valid)) {
    return(data.frame(
      network = character(), header = character(), value = character(),
      numeric_value = numeric(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    network = network,
    header = vapply(parts[valid], `[[`, character(1), 2L),
    value = vapply(parts[valid], `[[`, character(1), 3L),
    numeric_value = suppressWarnings(as.numeric(vapply(parts[valid], `[[`, character(1), 3L))),
    stringsAsFactors = FALSE
  )
}

#' Parse one hierarchical HotNet cluster result file
#'
#' Each non-comment line is retained as one reported cluster. Tab-separated
#' fields are the member node names; comment lines are parsed as result header
#' metadata.
#'
#' @param file Cluster TSV written by hierarchical HotNet.
#' @param network Network name used in output tables.
#' @return A list containing `headers`, `clusters`, and `membership` tables.
read_hhotnet_cluster_result <- function(file, network) {
  if (!file.exists(file)) stop("Missing h-HotNet cluster file: ", file, call. = FALSE)
  lines <- readLines(file, warn = FALSE)
  is_header <- grepl("^#", lines)
  headers <- parse_hhotnet_cluster_headers(lines[is_header], network)
  cluster_lines <- trimws(lines[!is_header])
  cluster_lines <- cluster_lines[nzchar(cluster_lines)]
  members <- lapply(cluster_lines, function(line) {
    values <- trimws(strsplit(line, "\\t", fixed = FALSE)[[1]])
    values[nzchar(values)]
  })
  if (!length(members)) stop("No clusters found in h-HotNet result: ", file, call. = FALSE)
  sizes <- lengths(members)
  clusters <- data.frame(
    network = network,
    cluster = seq_along(members),
    cluster_size = as.integer(sizes),
    is_singleton = sizes == 1L,
    is_multi_node = sizes > 1L,
    stringsAsFactors = FALSE
  )
  membership <- data.frame(
    network = network,
    cluster = rep(clusters$cluster, sizes),
    node = unlist(members, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(membership$node)) {
    duplicated_nodes <- unique(membership$node[duplicated(membership$node)])
    stop(
      "A h-HotNet node occurs in more than one reported cluster for ", network,
      ": ", paste(utils::head(duplicated_nodes, 10L), collapse = ", "), call. = FALSE
    )
  }
  list(headers = headers, clusters = clusters, membership = membership)
}

safe_hhotnet_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_hhotnet_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

hhotnet_header_field <- function(header) {
  field <- tolower(gsub("[^A-Za-z0-9]+", "_", trimws(header)))
  field <- gsub("(^_+|_+$)", "", field)
  paste0("header_", ifelse(nzchar(field), field, "unnamed"))
}

hhotnet_result_data <- function(config, network) {
  require_hhotnet_result_packages()
  paths <- analysis_paths(config)
  parsed <- read_hhotnet_cluster_result(hhotnet_cluster_result_file(config, network), network)
  edge_file <- network_data_file(config, paste0(network, "_edge_list_full.tsv"))
  node_file <- network_data_file(config, "nodes_all_info.tsv")
  if (!file.exists(edge_file) || !file.exists(node_file)) {
    stop("Missing h-HotNet network table(s) for ", network, call. = FALSE)
  }
  edges <- readr::read_tsv(edge_file, show_col_types = FALSE, progress = FALSE)
  required_edge_columns <- c("from", "to", "source")
  if (length(setdiff(required_edge_columns, names(edges)))) {
    stop("Network edge table lacks required columns: ", edge_file, call. = FALSE)
  }
  edges <- dplyr::mutate(
    edges,
    from = as.character(.data$from),
    to = as.character(.data$to),
    source = as.character(.data$source)
  )
  graph_nodes <- unique(c(edges$from, edges$to))
  nodes_all <- readr::read_tsv(node_file, show_col_types = FALSE, progress = FALSE)
  if (!"gene" %in% names(nodes_all)) stop("Node table lacks gene column: ", node_file, call. = FALSE)
  nodes_all$gene <- as.character(nodes_all$gene)
  if (anyDuplicated(nodes_all$gene)) stop("Node table contains duplicate gene identifiers: ", node_file, call. = FALSE)
  missing_graph_nodes <- setdiff(graph_nodes, nodes_all$gene)
  if (length(missing_graph_nodes)) {
    stop(
      "Node table lacks graph endpoint(s) for ", network, ": ",
      paste(utils::head(missing_graph_nodes, 10L), collapse = ", "), call. = FALSE
    )
  }
  unknown_cluster_nodes <- setdiff(parsed$membership$node, nodes_all$gene)
  if (length(unknown_cluster_nodes)) {
    stop(
      "Cluster result contains node(s) absent from nodes_all_info.tsv for ", network,
      ": ", paste(utils::head(unknown_cluster_nodes, 10L), collapse = ", "), call. = FALSE
    )
  }
  nodes <- nodes_all[match(graph_nodes, nodes_all$gene), , drop = FALSE]
  cluster_lookup <- dplyr::select(parsed$membership, "node", "cluster") |>
    dplyr::left_join(dplyr::select(parsed$clusters, "cluster", "cluster_size",
                                   "is_singleton", "is_multi_node"), by = "cluster")
  node_table <- dplyr::left_join(nodes, cluster_lookup, by = c("gene" = "node")) |>
    dplyr::mutate(
      cluster = dplyr::coalesce(as.integer(.data$cluster), 0L),
      cluster_size = dplyr::coalesce(as.integer(.data$cluster_size), 0L),
      is_singleton = dplyr::coalesce(.data$is_singleton, FALSE),
      is_multi_node = dplyr::coalesce(.data$is_multi_node, FALSE),
      in_hotnet_cluster = .data$is_multi_node,
      name = .data$gene,
      label = .data$gene
    ) |>
    dplyr::relocate("name", "label", "cluster", "cluster_size",
                    "is_singleton", "is_multi_node", "in_hotnet_cluster")
  member_node_attributes <- dplyr::select(
    node_table, -"cluster", -"cluster_size", -"is_singleton", -"is_multi_node",
    -"in_hotnet_cluster"
  )
  member_nodes <- dplyr::left_join(parsed$membership, member_node_attributes, by = c("node" = "name"))
  pvalue_column <- config$scores$pvalue_column
  cluster_summary <- member_nodes |>
    dplyr::group_by(.data$network, .data$cluster) |>
    dplyr::summarise(
      node_count = dplyr::n(),
      enhancer_count = sum(.data$node_type == "enhancer", na.rm = TRUE),
      promoter_gene_count = sum(.data$node_type == "promoter_gene", na.rm = TRUE),
      min_score = safe_hhotnet_min(.data$score),
      max_score = safe_hhotnet_max(.data$score),
      best_pvalue = safe_hhotnet_min(.data[[pvalue_column]]),
      genes = paste(.data$node[.data$node_type == "promoter_gene"], collapse = ";"),
      enhancers = paste(.data$node[.data$node_type == "enhancer"], collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::right_join(parsed$clusters, by = c("network", "cluster")) |>
    dplyr::mutate(
      node_count = dplyr::coalesce(as.integer(.data$node_count), 0L),
      enhancer_count = dplyr::coalesce(as.integer(.data$enhancer_count), 0L),
      promoter_gene_count = dplyr::coalesce(as.integer(.data$promoter_gene_count), 0L)
    ) |>
    dplyr::arrange(.data$cluster)
  membership_map <- stats::setNames(parsed$membership$cluster, parsed$membership$node)
  edge_cluster <- unname(membership_map[edges$from])
  target_cluster <- unname(membership_map[edges$to])
  same_cluster <- !is.na(edge_cluster) & edge_cluster == target_cluster
  cluster_edges <- edges[same_cluster, , drop = FALSE]
  cluster_edges <- dplyr::mutate(cluster_edges, network = network, cluster = as.integer(edge_cluster[same_cluster])) |>
    dplyr::relocate("network", "cluster")
  header_fields <- hhotnet_header_field(parsed$headers$header)
  header_fields <- make.unique(header_fields, sep = "_")
  header_values <- stats::setNames(as.list(parsed$headers$value), header_fields)
  result_summary <- data.frame(
    network = network,
    result_file = hhotnet_cluster_result_file(config, network),
    input_nodes = nrow(node_table),
    input_edges = nrow(edges),
    reported_clusters = nrow(parsed$clusters),
    multi_node_clusters = sum(parsed$clusters$is_multi_node),
    singleton_clusters = sum(parsed$clusters$is_singleton),
    clustered_nodes = nrow(parsed$membership),
    stringsAsFactors = FALSE
  )
  if (length(header_values)) result_summary <- cbind(result_summary, as.data.frame(header_values, stringsAsFactors = FALSE))
  list(
    headers = parsed$headers,
    clusters = cluster_summary,
    cluster_nodes = member_nodes,
    cluster_edges = cluster_edges,
    nodes = node_table,
    edges = edges,
    result_summary = result_summary
  )
}

summarize_hhotnet_permutation_outputs <- function(config, network) {
  require_hhotnet_result_packages()
  score_dir <- file.path(
    analysis_paths(config)[["hhotnet_intermediate"]],
    paste0(network, "_", hhotnet_score_name(config))
  )
  edge_files <- list.files(
    score_dir, pattern = "^hierarchy_edge_list_[0-9]+\\.tsv$", full.names = TRUE
  )
  if (!length(edge_files)) stop("No hierarchy edge files found for ", network, call. = FALSE)
  permutation <- as.integer(sub("^.*_([0-9]+)\\.tsv$", "\\1", edge_files))
  edge_files <- edge_files[order(permutation)]
  permutation <- sort(permutation)
  dplyr::bind_rows(lapply(seq_along(edge_files), function(i) {
    index_file <- file.path(score_dir, paste0("hierarchy_index_gene_", permutation[[i]], ".tsv"))
    if (!file.exists(index_file)) {
      stop("Missing hierarchy index file for ", network, " permutation ", permutation[[i]], call. = FALSE)
    }
    edges <- readr::read_tsv(
      edge_files[[i]], col_names = c("id1", "id2", "height"),
      col_types = readr::cols(id1 = readr::col_integer(), id2 = readr::col_integer(), height = readr::col_double()),
      show_col_types = FALSE, progress = FALSE
    )
    index <- readr::read_tsv(
      index_file, col_names = c("id", "gene"),
      col_types = readr::cols(id = readr::col_integer(), gene = readr::col_character()),
      show_col_types = FALSE, progress = FALSE
    )
    valid_index <- identical(index$id, seq_len(nrow(index))) && !anyDuplicated(index$gene)
    # Hierarchy edge lists contain both leaf IDs from the gene index and new
    # internal dendrogram IDs. Internal IDs are expected to exceed nrow(index).
    positive_hierarchy_ids <- all(is.finite(edges$id1) & edges$id1 > 0L) &&
      all(is.finite(edges$id2) & edges$id2 > 0L)
    edge_endpoints_in_gene_index <- all(edges$id1 %in% index$id) && all(edges$id2 %in% index$id)
    data.frame(
      network = network,
      permutation = permutation[[i]],
      hierarchy_kind = if (permutation[[i]] == 0L) "observed" else "permuted",
      index_nodes = nrow(index), hierarchy_edges = nrow(edges),
      min_height = safe_hhotnet_min(edges$height), max_height = safe_hhotnet_max(edges$height),
      valid_index = valid_index, positive_hierarchy_ids = positive_hierarchy_ids,
      edge_endpoints_in_gene_index = edge_endpoints_in_gene_index,
      edge_file = edge_files[[i]], index_file = index_file,
      stringsAsFactors = FALSE
    )
  }))
}

write_hhotnet_cluster_tables <- function(result, config, network) {
  cluster_root <- file.path(analysis_paths(config)[["hhotnet_summary"]], "clusters", network)
  dir.create(cluster_root, recursive = TRUE, showWarnings = FALSE)
  selected <- dplyr::filter(result$clusters, .data$is_multi_node)
  manifest <- lapply(seq_len(nrow(selected)), function(i) {
    cluster_id <- selected$cluster[[i]]
    stem <- sprintf("cluster_%03d", cluster_id)
    node_file <- file.path(cluster_root, paste0(stem, "_nodes.tsv"))
    edge_file <- file.path(cluster_root, paste0(stem, "_edges.tsv"))
    readr::write_tsv(
      result$cluster_nodes[result$cluster_nodes$cluster == cluster_id, , drop = FALSE],
      node_file
    )
    readr::write_tsv(
      result$cluster_edges[result$cluster_edges$cluster == cluster_id, , drop = FALSE],
      edge_file
    )
    data.frame(network = network, cluster = cluster_id, node_file = node_file,
               edge_file = edge_file, stringsAsFactors = FALSE)
  })
  manifest <- if (length(manifest)) dplyr::bind_rows(manifest) else data.frame(
    network = character(), cluster = integer(), node_file = character(), edge_file = character()
  )
  manifest_file <- file.path(cluster_root, "cluster_table_manifest.tsv")
  readr::write_tsv(manifest, manifest_file)
  manifest_file
}

#' Summarize completed hierarchical HotNet results
#'
#' Parses result headers and every reported cluster, joins cluster members to
#' the exact input node table, and writes complete cluster-specific induced edge
#' tables. No h-HotNet calculation is rerun.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for the configured default, one network name, or a
#'   vector of network names.
#' @return A network-level result summary data frame.
#' @export
summarize_hhotnet_results <- function(config, networks = NULL) {
  require_hhotnet_result_packages()
  requested <- resolve_networks(config, networks)
  paths <- analysis_paths(config)
  dir.create(paths[["hhotnet_summary"]], recursive = TRUE, showWarnings = FALSE)
  results <- lapply(requested, function(network) hhotnet_result_data(config, network))
  names(results) <- requested
  permutations <- lapply(requested, function(network) summarize_hhotnet_permutation_outputs(config, network))
  names(permutations) <- requested
  for (network in requested) {
    result <- results[[network]]
    readr::write_tsv(result$headers, network_summary_file(config, paste0(network, "_headers.tsv")))
    readr::write_tsv(result$clusters, network_summary_file(config, paste0(network, "_cluster_summary.tsv")))
    readr::write_tsv(result$cluster_nodes, network_summary_file(config, paste0(network, "_cluster_nodes.tsv")))
    readr::write_tsv(result$cluster_edges, network_summary_file(config, paste0(network, "_cluster_edges.tsv")))
    readr::write_tsv(result$result_summary, network_summary_file(config, paste0(network, "_result_summary.tsv")))
    readr::write_tsv(permutations[[network]], network_summary_file(config, paste0(network, "_permutation_summary.tsv")))
    write_hhotnet_cluster_tables(result, config, network)
  }
  combined <- list(
    headers = dplyr::bind_rows(lapply(results, `[[`, "headers")),
    clusters = dplyr::bind_rows(lapply(results, `[[`, "clusters")),
    cluster_nodes = dplyr::bind_rows(lapply(results, `[[`, "cluster_nodes")),
    cluster_edges = dplyr::bind_rows(lapply(results, `[[`, "cluster_edges")),
    result_summary = dplyr::bind_rows(lapply(results, `[[`, "result_summary")),
    permutations = dplyr::bind_rows(permutations)
  )
  readr::write_tsv(combined$headers, network_summary_file(config, "all_network_headers.tsv"))
  readr::write_tsv(combined$clusters, network_summary_file(config, "all_network_cluster_summary.tsv"))
  readr::write_tsv(combined$cluster_nodes, network_summary_file(config, "all_network_cluster_nodes.tsv"))
  readr::write_tsv(combined$cluster_edges, network_summary_file(config, "all_network_cluster_edges.tsv"))
  readr::write_tsv(combined$result_summary, network_summary_file(config, "all_network_result_summary.tsv"))
  readr::write_tsv(combined$permutations, network_summary_file(config, "all_network_permutation_summary.tsv"))
  report <- c(
    "# h-HotNet result summary",
    "",
    paste0("Analysis: `", config$analysis$name, "`."),
    paste0("Configuration: `", config$config_file, "`."),
    paste0("Score table: `", hhotnet_score_name(config), "`; score permutations: ",
           config$hhotnet$num_permutations, "; lower cluster-size bound: ",
           config$hhotnet$lower_size_bound, "."),
    "",
    "| Network | Input nodes | Input edges | Reported clusters | Multi-node clusters | Singleton clusters |",
    "| --- | ---: | ---: | ---: | ---: | ---: |",
    vapply(seq_len(nrow(combined$result_summary)), function(i) {
      row <- combined$result_summary[i, ]
      paste0("| ", row$network, " | ", row$input_nodes, " | ", row$input_edges,
             " | ", row$reported_clusters, " | ", row$multi_node_clusters,
             " | ", row$singleton_clusters, " |")
    }, character(1)),
    "",
    "Each non-comment line in the h-HotNet cluster file is represented as one cluster.",
    "The permutation summary tabulates the observed hierarchy (permutation 0) and every completed score permutation.",
    "Cluster-specific node and induced-edge TSV files are written for multi-node clusters."
  )
  writeLines(report, network_summary_file(config, "hhotnet_result_report.md"))
  log_analysis_event(config, "hhotnet_results_summarized", paste(requested, collapse = ","))
  combined$result_summary
}

scale_hhotnet_node_size <- function(score, size_range = c(20, 60)) {
  score <- suppressWarnings(as.numeric(score))
  finite <- is.finite(score)
  output <- rep(mean(size_range), length(score))
  if (!any(finite)) return(output)
  limits <- range(score[finite])
  if (diff(limits) == 0) return(output)
  output[finite] <- size_range[[1]] +
    (score[finite] - limits[[1]]) / diff(limits) * diff(size_range)
  output
}

style_hhotnet_nodes <- function(nodes, size_range = c(20, 60)) {
  dplyr::mutate(
    nodes,
    node_class = dplyr::if_else(.data$node_type == "promoter_gene", "gene", "enhancer"),
    background_source = dplyr::case_when(
      grepl("JEME|HiC", .data$node_sources) ~ "JEME/HiC",
      grepl("ARACNe", .data$node_sources) ~ "ARACNe only",
      TRUE ~ "other"
    ),
    shape = dplyr::if_else(.data$node_class == "gene", "DIAMOND", "ELLIPSE"),
    cy_shape = dplyr::if_else(.data$node_class == "gene", "diamond", "ellipse"),
    fill_color = dplyr::case_when(
      .data$in_hotnet_cluster & .data$node_class == "gene" ~ "#B2182B",
      .data$in_hotnet_cluster & .data$node_class == "enhancer" ~ "#EF8A8A",
      .data$background_source == "JEME/HiC" ~ "#ADD8E6",
      TRUE ~ "#D3D3D3"
    ),
    node_size = scale_hhotnet_node_size(.data$score, size_range)
  )
}

style_hhotnet_edges <- function(edges) {
  dplyr::mutate(
    edges,
    interaction = .data$source,
    edge_class = dplyr::case_when(
      .data$source == "JEME" ~ "JEME",
      grepl("^HiC", .data$source) ~ "HiC",
      .data$source == "ARACNe" ~ "ARACNe",
      TRUE ~ "other"
    ),
    edge_color = dplyr::case_when(
      .data$edge_class == "JEME" ~ "#4C78A8",
      .data$edge_class == "HiC" ~ "#F2A541",
      .data$edge_class == "ARACNe" ~ "#9E9E9E",
      TRUE ~ "#D3D3D3"
    ),
    edge_width = dplyr::if_else(.data$edge_class %in% c("JEME", "HiC"), 1.5, 1),
    edge_line_style = dplyr::if_else(.data$edge_class == "ARACNe", "DOT", "SOLID")
  )
}

filter_hhotnet_component_tables <- function(nodes, edges) {
  graph <- igraph::graph_from_data_frame(
    dplyr::select(edges, "from", "to"), directed = FALSE,
    vertices = data.frame(name = nodes$name, stringsAsFactors = FALSE)
  )
  components <- igraph::components(graph)
  memberships <- stats::setNames(components$membership, igraph::V(graph)$name)
  cluster_components <- unique(memberships[nodes$name[nodes$in_hotnet_cluster]])
  cluster_components <- cluster_components[!is.na(cluster_components)]
  if (!length(cluster_components)) {
    return(list(nodes = nodes[0, , drop = FALSE], edges = edges[0, , drop = FALSE]))
  }
  keep_names <- names(memberships)[memberships %in% cluster_components]
  list(
    nodes = dplyr::filter(nodes, .data$name %in% keep_names),
    edges = dplyr::filter(edges, .data$from %in% keep_names, .data$to %in% keep_names)
  )
}

hhotnet_igraph_from_tables <- function(nodes, edges) {
  edge_columns <- c("from", "to", setdiff(names(edges), c("from", "to")))
  igraph::graph_from_data_frame(edges[, edge_columns, drop = FALSE], directed = FALSE, vertices = nodes)
}

row_as_json_list <- function(data, index) {
  values <- as.list(data[index, , drop = FALSE])
  lapply(values, function(value) unname(value[[1]]))
}

cx2_type <- function(x) {
  if (is.logical(x)) return("boolean")
  if (is.integer(x)) return("integer")
  if (is.numeric(x)) return("double")
  "string"
}

cx2_attribute_declarations <- function(data) {
  stats::setNames(lapply(data, function(column) list(d = cx2_type(column))), names(data))
}

write_hhotnet_cytoscape_json <- function(nodes, edges, file, analysis, network) {
  node_elements <- lapply(seq_len(nrow(nodes)), function(i) {
    data <- row_as_json_list(nodes, i)
    data$id <- data$name
    list(data = data)
  })
  edge_elements <- lapply(seq_len(nrow(edges)), function(i) {
    data <- row_as_json_list(edges, i)
    data$id <- paste0("edge_", i)
    data$network_source <- data$source
    data$source <- data$from
    data$target <- data$to
    list(data = data)
  })
  cyjs <- list(
    format_version = "1.0",
    generated_by = "hotnetR2::export_hhotnet_graphs",
    analysis = analysis,
    network = network,
    elements = list(nodes = node_elements, edges = edge_elements),
    style = list(
      list(selector = "node", style = list(
        label = "data(label)", shape = "data(cy_shape)", width = "data(node_size)",
        height = "data(node_size)", `background-color` = "data(fill_color)", `border-width` = 0
      )),
      list(selector = "edge", style = list(
        width = "data(edge_width)", `line-color` = "data(edge_color)",
        `line-style` = "data(edge_line_style)", `curve-style` = "bezier"
      ))
    )
  )
  jsonlite::write_json(cyjs, file, pretty = TRUE, auto_unbox = TRUE, na = "null")
}

write_hhotnet_cx2 <- function(nodes, edges, file, analysis, network, component_filter) {
  network_name <- tools::file_path_sans_ext(basename(file))
  node_ids <- stats::setNames(seq_len(nrow(nodes)) - 1L, nodes$name)
  if (anyDuplicated(names(node_ids))) stop("CX2 node names must be unique", call. = FALSE)
  cx_nodes <- lapply(seq_len(nrow(nodes)), function(i) {
    list(id = unname(node_ids[[nodes$name[[i]]]]), v = row_as_json_list(nodes, i))
  })
  cx_edges <- lapply(seq_len(nrow(edges)), function(i) {
    list(
      id = i - 1L,
      s = unname(node_ids[[edges$from[[i]]]]),
      t = unname(node_ids[[edges$to[[i]]]]),
      v = row_as_json_list(edges, i)
    )
  })
  node_bypasses <- lapply(seq_len(nrow(nodes)), function(i) {
    list(id = unname(node_ids[[nodes$name[[i]]]]), v = list(
      NODE_SHAPE = nodes$cy_shape[[i]], NODE_BACKGROUND_COLOR = nodes$fill_color[[i]],
      NODE_WIDTH = nodes$node_size[[i]], NODE_HEIGHT = nodes$node_size[[i]]
    ))
  })
  edge_bypasses <- lapply(seq_len(nrow(edges)), function(i) {
    list(id = i - 1L, v = list(
      EDGE_LINE_COLOR = edges$edge_color[[i]], EDGE_WIDTH = edges$edge_width[[i]],
      EDGE_LINE_STYLE = edges$edge_line_style[[i]]
    ))
  })
  cx2 <- list(
    list(CXVersion = "2.0", hasFragments = FALSE),
    list(attributeDeclarations = list(list(
      networkAttributes = list(
        name = list(d = "string"), analysis = list(d = "string"), network = list(d = "string"),
        component_filter = list(d = "boolean")
      ),
      nodes = cx2_attribute_declarations(nodes),
      edges = cx2_attribute_declarations(edges)
    ))),
    list(networkAttributes = list(list(
      name = network_name, analysis = analysis, network = network,
      component_filter = isTRUE(component_filter)
    ))),
    list(nodes = cx_nodes),
    list(edges = cx_edges),
    list(visualProperties = list(list(
      default = list(
        network = list(NETWORK_BACKGROUND_COLOR = "#FFFFFF"),
        node = list(NODE_SHAPE = "ELLIPSE", NODE_BACKGROUND_COLOR = "#D3D3D3", NODE_WIDTH = 30, NODE_HEIGHT = 30),
        edge = list(EDGE_WIDTH = 1, EDGE_LINE_COLOR = "#D3D3D3")
      ),
      nodeMapping = list(NODE_LABEL = list(
        type = "PASSTHROUGH", definition = list(attribute = "label", type = "string")
      ))
    ))),
    list(nodeBypasses = node_bypasses),
    list(edgeBypasses = edge_bypasses),
    list(visualEditorProperties = list(list(nodeSizeLocked = TRUE))),
    list(status = list(list(success = TRUE, error = "")))
  )
  jsonlite::write_json(cx2, file, pretty = TRUE, auto_unbox = TRUE, na = "null")
}

hhotnet_export_prefix <- function(config, network, filter_components) {
  suffix <- if (isTRUE(filter_components)) "hotnet_components" else "all"
  paste(config$analysis$name, network, suffix, sep = "_")
}

#' Export completed h-HotNet input graphs and attributes
#'
#' Produces a complete node TSV, complete edge TSV, GraphML, Cytoscape.js JSON,
#' and CX2. CX2 preserves every column of the corresponding exported edge TSV.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for the configured default, one network name, or a
#'   vector of network names.
#' @param filter_components If `TRUE`, retain connected components containing a
#'   reported multi-node h-HotNet cluster.
#' @param size_range Two increasing display sizes for node-score scaling.
#' @return A data frame listing written exports.
#' @export
export_hhotnet_graphs <- function(config, networks = NULL, filter_components = TRUE,
                                  size_range = c(20, 60)) {
  require_hhotnet_result_packages(include_json = TRUE)
  requested <- resolve_networks(config, networks)
  if (length(size_range) != 2L || any(!is.finite(size_range)) || size_range[[1]] >= size_range[[2]]) {
    stop("size_range must contain two increasing finite values", call. = FALSE)
  }
  paths <- analysis_paths(config)
  dir.create(paths[["hhotnet_cytoscape"]], recursive = TRUE, showWarnings = FALSE)
  manifest <- lapply(requested, function(network) {
    result <- hhotnet_result_data(config, network)
    nodes <- style_hhotnet_nodes(result$nodes, size_range)
    edges <- style_hhotnet_edges(result$edges)
    selected <- if (isTRUE(filter_components)) filter_hhotnet_component_tables(nodes, edges) else list(nodes = nodes, edges = edges)
    if (!nrow(selected$nodes) || !nrow(selected$edges)) {
      stop("No graph content remains for export after filtering ", network, call. = FALSE)
    }
    prefix <- hhotnet_export_prefix(config, network, filter_components)
    node_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, "_nodes.tsv"))
    edge_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, "_edges.tsv"))
    graphml_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, ".graphml"))
    cyjs_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, ".cyjs.json"))
    cx2_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, ".cx2"))
    readr::write_tsv(selected$nodes, node_file, na = "NA")
    readr::write_tsv(selected$edges, edge_file, na = "NA")
    graph <- hhotnet_igraph_from_tables(selected$nodes, selected$edges)
    igraph::write_graph(graph, graphml_file, format = "graphml")
    write_hhotnet_cytoscape_json(selected$nodes, selected$edges, cyjs_file, config$analysis$name, network)
    write_hhotnet_cx2(selected$nodes, selected$edges, cx2_file, config$analysis$name, network, filter_components)
    data.frame(
      network = network, filter_components = isTRUE(filter_components), nodes = nrow(selected$nodes), edges = nrow(selected$edges),
      node_file = node_file, edge_file = edge_file, graphml_file = graphml_file,
      cytoscape_json_file = cyjs_file, cx2_file = cx2_file, stringsAsFactors = FALSE
    )
  })
  manifest <- dplyr::bind_rows(manifest)
  readr::write_tsv(manifest, file.path(paths[["hhotnet_cytoscape"]], "export_manifest.tsv"))
  log_analysis_event(config, "hhotnet_graphs_exported", paste(requested, collapse = ","))
  manifest
}

cx2_aspect <- function(cx2, name) {
  index <- which(vapply(cx2, function(aspect) name %in% names(aspect), logical(1)))
  if (!length(index)) return(NULL)
  cx2[[index[[1]]]][[name]]
}

#' Validate Phase 8 graph exports
#'
#' Checks TSV endpoint integrity, GraphML readability, Cytoscape.js element
#' counts, and CX2 structural IDs and edge-attribute declarations.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for the configured default, one network name, or a
#'   vector of network names.
#' @param filter_components Match the setting used when graph exports were written.
#' @return A validation table, written below `hHotnet/summary`.
#' @export
validate_hhotnet_exports <- function(config, networks = NULL, filter_components = TRUE) {
  require_hhotnet_result_packages(include_json = TRUE)
  requested <- resolve_networks(config, networks)
  paths <- analysis_paths(config)
  report <- lapply(requested, function(network) {
    prefix <- hhotnet_export_prefix(config, network, filter_components)
    node_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, "_nodes.tsv"))
    edge_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, "_edges.tsv"))
    graphml_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, ".graphml"))
    cyjs_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, ".cyjs.json"))
    cx2_file <- file.path(paths[["hhotnet_cytoscape"]], paste0(prefix, ".cx2"))
    files_exist <- all(file.exists(c(node_file, edge_file, graphml_file, cyjs_file, cx2_file)))
    if (!files_exist) {
      return(data.frame(network = network, files_exist = FALSE, unique_node_names = FALSE,
                        valid_endpoints = FALSE, graphml_readable = FALSE,
                        cyjs_counts_match = FALSE, cx2_counts_match = FALSE,
                        cx2_valid_endpoints = FALSE, cx2_all_edge_attributes = FALSE,
                        valid = FALSE, stringsAsFactors = FALSE))
    }
    text_columns <- readr::cols(.default = readr::col_character())
    nodes <- readr::read_tsv(node_file, col_types = text_columns, show_col_types = FALSE, progress = FALSE)
    edges <- readr::read_tsv(edge_file, col_types = text_columns, show_col_types = FALSE, progress = FALSE)
    unique_node_names <- "name" %in% names(nodes) && !anyDuplicated(nodes$name) && all(nzchar(nodes$name))
    valid_endpoints <- all(c("from", "to") %in% names(edges)) && all(edges$from %in% nodes$name) && all(edges$to %in% nodes$name)
    graphml_readable <- !inherits(try(igraph::read_graph(graphml_file, format = "graphml"), silent = TRUE), "try-error")
    cyjs <- jsonlite::read_json(cyjs_file, simplifyVector = FALSE)
    cyjs_counts_match <- !is.null(cyjs$elements) && length(cyjs$elements$nodes) == nrow(nodes) && length(cyjs$elements$edges) == nrow(edges)
    cx2 <- jsonlite::read_json(cx2_file, simplifyVector = FALSE)
    cx_nodes <- cx2_aspect(cx2, "nodes")
    cx_edges <- cx2_aspect(cx2, "edges")
    declarations <- cx2_aspect(cx2, "attributeDeclarations")
    cx2_counts_match <- length(cx_nodes) == nrow(nodes) && length(cx_edges) == nrow(edges)
    node_ids <- if (length(cx_nodes)) vapply(cx_nodes, function(x) as.integer(x$id), integer(1)) else integer()
    edge_ids <- if (length(cx_edges)) vapply(cx_edges, function(x) as.integer(x$id), integer(1)) else integer()
    endpoint_ids <- if (length(cx_edges)) unlist(lapply(cx_edges, function(x) c(as.integer(x$s), as.integer(x$t)))) else integer()
    cx2_valid_endpoints <- !anyDuplicated(node_ids) && !anyDuplicated(edge_ids) && all(endpoint_ids %in% node_ids)
    declared_edges <- if (!is.null(declarations) && length(declarations) && !is.null(declarations[[1]]$edges)) names(declarations[[1]]$edges) else character()
    cx2_edge_values_complete <- if (length(cx_edges)) {
      all(vapply(cx_edges, function(edge) all(names(edges) %in% names(edge$v)), logical(1)))
    } else TRUE
    cx2_all_edge_attributes <- all(names(edges) %in% declared_edges) && cx2_edge_values_complete
    valid <- files_exist && unique_node_names && valid_endpoints && graphml_readable && cyjs_counts_match && cx2_counts_match && cx2_valid_endpoints && cx2_all_edge_attributes
    data.frame(network = network, files_exist = files_exist, unique_node_names = unique_node_names,
               valid_endpoints = valid_endpoints, graphml_readable = graphml_readable,
               cyjs_counts_match = cyjs_counts_match, cx2_counts_match = cx2_counts_match,
               cx2_valid_endpoints = cx2_valid_endpoints, cx2_all_edge_attributes = cx2_all_edge_attributes,
               valid = valid, stringsAsFactors = FALSE)
  })
  report <- dplyr::bind_rows(report)
  readr::write_tsv(report, network_summary_file(config, "hhotnet_export_validation.tsv"))
  if (any(!report$valid)) stop("One or more h-HotNet exports failed validation", call. = FALSE)
  report
}
