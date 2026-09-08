#' Add node attributes to an igraph object
#'
#' This function adds attributes from a data frame to the vertices of an `igraph` object.
#' Optionally, it can remove vertices not present in the data frame and ensure only nodes
#' present in the graph are used for assignment.
#'
#' @param ig An `igraph` object. Vertices must have a `name` attribute.
#' @param nodes A data frame or tibble containing node attributes. One column must match
#'   the vertex names in the graph (default is `"gene"`).
#' @param node_name Character string. Column name in `nodes` used to match graph vertices (default = `"gene"`).
#' @param rm_miss Logical. If `TRUE`, remove vertices in the graph not present in `nodes` (default = `TRUE`).
#' @param reverse_check Logical. If `TRUE`, filter `nodes` to keep only rows present in the graph (default = `TRUE`).
#'
#' @return An `igraph` object with updated vertex attributes.
#' @importFrom igraph delete_vertices V set_vertex_attr vertex_attr
#' @importFrom dplyr filter
#' @export
add_node_attributes <- function(ig, nodes, node_name = "gene",
                                rm_miss = TRUE, reverse_check = TRUE) {
  stopifnot(node_name %in% names(nodes))

  node_map <- nodes[[node_name]]

  # Remove vertices not found in nodes
  if (rm_miss) {
    ig <- igraph::delete_vertices(ig, setdiff(igraph::V(ig)$name, node_map))
  }

  # Optionally check if all nodes dataframe are present in ig
  if (reverse_check) {
    nodes <- nodes |> dplyr::filter(.data[[node_name]] %in% igraph::V(ig)$name)
  }

  # Assign attributes
  for (attr in setdiff(names(nodes), node_name)) {
    vals <- stats::setNames(nodes[[attr]], nodes[[node_name]])
    ig <- igraph::set_vertex_attr(ig, attr, index = names(vals), value = vals)
  }

  ig
}

#' Identify the vertex with minimum score in each subgraph/component
#'
#' This function computes the minimum value of a numeric vertex attribute (e.g., "score")
#' within each connected component (subgraph) of an igraph object and returns a tidy tibble
#' summarizing the result.
#'
#' @param ig An `igraph` object. Can be directed or undirected.
#' @param attr Character string. Name of the numeric vertex attribute to evaluate (default = "score").
#'
#' @return A tibble with one row per connected component, containing:
#' \describe{
#'   \item{component}{Component ID (integer).}
#'   \item{vertex}{Vertex name with the minimum score in that component (or numeric ID if names are not set).}
#'   \item{min_score}{Minimum score value within that component.}
#' }
#' @importFrom igraph components vertex_attr V
#' @importFrom tibble tibble
#' @export
min_score_per_component <- function(ig, attr = "score") {
  comp <- igraph::components(ig)
  scores <- igraph::vertex_attr(ig, attr)

  tibble::tibble(
    component = seq_len(comp$no),
    vertex = sapply(seq_len(comp$no), function(cid) {
      nodes <- which(comp$membership == cid)
      v <- nodes[which.min(scores[nodes])]
      if (!is.null(igraph::V(ig)$name)) igraph::V(ig)$name[v] else v
    }),
    min_score = sapply(seq_len(comp$no), function(cid) {
      nodes <- which(comp$membership == cid)
      min(scores[nodes], na.rm = TRUE)
    })
  )
}

#' Count vertices and edges in each connected component of a graph
#'
#' This function decomposes an igraph object into its connected components (subgraphs)
#' and returns a summary table containing the number of vertices and edges in each component.
#'
#' @param g An `igraph` object. Can be directed or undirected.
#'
#' @return A data frame with one row per connected component containing:
#' \describe{
#'   \item{component}{Component ID (integer).}
#'   \item{nodes}{Number of vertices in the component.}
#'   \item{edges}{Number of edges in the component.}
#' }
#' @importFrom igraph decompose gorder gsize
#' @export
count_edges_per_component <- function(g) {
  # Decompose into connected components
  comps <- igraph::decompose(g)

  # Build summary table
  summary_df <- data.frame(
    component = seq_along(comps),
    nodes     = sapply(comps, igraph::gorder),  # number of vertices
    edges     = sapply(comps, igraph::gsize)    # number of edges
  )

  return(summary_df)
}

#' Summarize subgraph/component information in an igraph object
#'
#' This function decomposes an igraph object into connected components (subgraphs)
#' and returns a tibble summarizing each component, including:
#' the vertex with the maximum score, the maximum score itself, and the number
#' of nodes and edges in the component.
#'
#' @param ig An `igraph` object. Vertices must have a numeric attribute specified by `score_attr`.
#' @param score_attr Character string. Name of the numeric vertex attribute used to determine max score per component (default = `"score"`).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{component}{Component ID (integer).}
#'   \item{vertex}{Vertex name with the maximum score in that component.}
#'   \item{max_score}{Maximum score in the component.}
#'   \item{n_nodes}{Number of nodes in the component.}
#'   \item{n_edges}{Number of edges in the component.}
#' }
#' @importFrom igraph components V vertex_attr induced_subgraph vcount ecount
#' @importFrom tibble tibble
#' @importFrom dplyr group_by slice_max rename ungroup left_join select
#' @importFrom purrr map_dfr
#' @export
subgraph_info <- function(ig, score_attr = "score") {
  comps <- igraph::components(ig)
  igraph::V(ig)$component <- comps$membership

  # Vertex data
  v_df <- tibble::tibble(
    vertex = names(igraph::V(ig)),
    component = igraph::V(ig)$component,
    score = as.numeric(igraph::vertex_attr(ig, score_attr))
  )

  # Max score per component
  max_scores <- v_df |>
    dplyr::group_by(component) |>
    dplyr::slice_max(order_by = score, n = 1, with_ties = FALSE) |>
    dplyr::rename(max_score = score) |>
    dplyr::ungroup()

  # Node/edge counts per component
  comp_info <- purrr::map_dfr(unique(comps$membership), function(comp) {
    sg <- igraph::induced_subgraph(ig, vids = which(comps$membership == comp))
    tibble::tibble(
      component = comp,
      n_nodes   = igraph::vcount(sg),
      n_edges   = igraph::ecount(sg)
    )
  })

  # Join results
  result <- max_scores |>
    dplyr::left_join(comp_info, by = "component") |>
    dplyr::select(component, vertex, max_score, n_nodes, n_edges)

  return(result)
}

#' Summarize subgraph/component information in an igraph object (using FDR)
#'
#' @param ig An `igraph` object.
#' @param score_attr Character string.
#' @return A tibble with columns.
#' @importFrom igraph components V vertex_attr induced_subgraph vcount ecount
#' @importFrom tibble tibble
#' @importFrom dplyr group_by slice_min rename ungroup left_join select
#' @importFrom purrr map_dfr
#' @export
subgraph_info_FDR <- function(ig, score_attr = "score") {
  comps <- igraph::components(ig)
  igraph::V(ig)$component <- comps$membership

  # Vertex data
  v_df <- tibble::tibble(
    vertex = names(igraph::V(ig)),
    component = igraph::V(ig)$component,
    score = as.numeric(igraph::vertex_attr(ig, score_attr))
  )

  # Min score per component
  min_scores <- v_df |>
    dplyr::group_by(component) |>
    dplyr::slice_min(order_by = score, n = 1, with_ties = FALSE) |>
    dplyr::rename(min_score = score) |>
    dplyr::ungroup()

  # Node/edge counts per component
  comp_info <- purrr::map_dfr(unique(comps$membership), function(comp) {
    sg <- igraph::induced_subgraph(ig, vids = which(comps$membership == comp))
    tibble::tibble(
      component = comp,
      n_nodes   = igraph::vcount(sg),
      n_edges   = igraph::ecount(sg)
    )
  })

  # Join results
  result <- min_scores |>
    dplyr::left_join(comp_info, by = "component") |>
    dplyr::select(component, vertex, min_score, n_nodes, n_edges)

  return(result)
}

#' Subset an igraph object by component maximum vertex score
#'
#' This function splits an igraph object into connected components and retains only
#' components whose maximum vertex score exceeds a specified threshold. The resulting
#' subgraph includes only vertices in these high-scoring components.
#'
#' @param ig An `igraph` object. Vertices must have a numeric attribute specified by `score_attr`.
#' @param score_attr Character string. Name of the numeric vertex attribute to evaluate per component (default = `"score"`).
#' @param max_score_threshold Numeric. Minimum maximum score required for a component to be retained (default = 4).
#'
#' @return An igraph object representing the induced subgraph containing only components
#'   whose maximum vertex score exceeds `max_score_threshold`.
#' @importFrom igraph components V vertex_attr induced_subgraph
#' @export
subset_graph_by_score <- function(ig, score_attr = "score", max_score_threshold = 4) {
  # Compute components
  comps <- igraph::components(ig)
  igraph::V(ig)$component <- comps$membership

  # Vertex scores
  scores <- igraph::vertex_attr(ig, score_attr)

  # Max score per component
  max_scores <- sapply(unique(comps$membership), function(cid) {
    nodes <- which(comps$membership == cid)
    max(scores[nodes], na.rm = TRUE)
  })
  names(max_scores) <- unique(comps$membership)

  # Components to keep
  keep_components <- as.integer(names(max_scores)[max_scores >= max_score_threshold])

  # Vertices to keep
  keep_vertices <- igraph::V(ig)[igraph::V(ig)$component %in% keep_components]

  # Return induced subgraph
  igraph::induced_subgraph(ig, vids = keep_vertices)
}

#' Convert an igraph object to HotNet-compatible format
#'
#' This function takes an igraph object and prepares data structures suitable for HotNet analysis.
#' It creates a mapping between vertex names and integer IDs, extracts vertex attributes,
#' and converts edges to use integer IDs instead of vertex names.
#'
#' @param ig An `igraph` object.
#'
#' @return A list containing ig, scores, index, edge, and edge_full.
#' @importFrom igraph V as_data_frame
#' @importFrom tibble tibble
#' @importFrom dplyr mutate relocate rename inner_join select everything
#' @export
igraph_to_hotnet <- function(ig){
  # create index-to-vertex list
  idx <- tibble::tibble(vertex = igraph::V(ig)$name) |>
    dplyr::mutate(id = dplyr::row_number()) |>
    dplyr::relocate(id, .before = 1)

  # Vertex scores
  scores <- igraph::as_data_frame(ig, what = "vertices")

  # Edge list with integer indices
  edge  <- igraph::as_data_frame(ig, what = "edges") |>
    dplyr::inner_join(idx |> dplyr::rename(from = vertex, id1 = id), by = "from") |>
    dplyr::inner_join(idx |> dplyr::rename(to = vertex, id2 = id), by = "to") |>
    dplyr::select(id1, id2, dplyr::everything())

  return(list(ig = ig, scores = scores, index = idx,
      edge = edge |> dplyr::select(id1, id2), edge_full = edge))
}

#' Plot each connected component of an igraph object to a multi-page PDF
#'
#' @param ig An `igraph` object.
#' @param score_attr Character string.
#' @param pdf_file Character string.
#'
#' @return No return value. A PDF file is created.
#' @importFrom igraph vertex_attr_names components induced_subgraph V vcount
#' @export
plot_subgraphs_to_pdf <- function(ig, score_attr = "score", pdf_file = "subgraphs.pdf") {
  # Ensure the score attribute exists
  if (!score_attr %in% igraph::vertex_attr_names(ig)) {
    stop(paste("Vertex attribute", score_attr, "not found in the graph"))
  }

  # Get components
  comps <- igraph::components(ig)

  # Split vertices by component
  subgraphs <- lapply(unique(comps$membership), function(cid) {
    igraph::induced_subgraph(ig, vids = which(comps$membership == cid))
  })

  # Open PDF device
  grDevices::pdf(pdf_file)

  # Plot each subgraph
  for (sg in subgraphs) {
    # Node size proportional to score (optional scaling)
    node_score <- igraph::vertex_attr(sg, score_attr)
    node_size <- 5 + 10 * (node_score / max(node_score, na.rm = TRUE))  # scale to 5-15

    graphics::plot(
      sg,
      vertex.size = node_size,
      vertex.label = igraph::V(sg)$name,
      main = paste("Component with", igraph::vcount(sg), "nodes"),
      vertex.color = "skyblue",
      edge.color = "gray50",
      layout = igraph::layout_with_fr
    )
  }

  # Close PDF device
  grDevices::dev.off()

  message("Subgraphs saved to ", pdf_file)
}

#' Print summary information of an igraph network
#'
#' @param ig An igraph object.
#' @param score_attr Character string.
#' @param top_n Number of top/bottom components to display.
#'
#' @return Invisibly returns a tibble with subgraph/component information.
#' @importFrom igraph vcount ecount
#' @importFrom dplyr arrange desc
#' @export
print_network_info <- function(ig, score_attr = "score", top_n = 6) {
  cat("Number of vertices:", igraph::vcount(ig), "\n")
  cat("Number of edges:", igraph::ecount(ig), "\n\n")

  info <- subgraph_info(ig, score_attr = score_attr)
  cat("Number of connected components:", nrow(info), "\n\n")

  cat("Top", top_n, "components by number of edges:\n")
  print(dplyr::arrange(info, dplyr::desc(n_edges)) |> utils::head(top_n))

  cat("\nBottom", top_n, "components by number of edges:\n")
  print(dplyr::arrange(info, n_edges) |> utils::head(top_n))

  cat("\nTop", top_n, "components by max score:\n")
  print(dplyr::arrange(info, dplyr::desc(max_score)) |> utils::head(top_n))

  cat("\nBottom", top_n, "components by max score:\n")
  print(dplyr::arrange(info, max_score) |> utils::head(top_n))

  invisible(info)
}

#' Add ARACNe edges to an existing igraph network
#'
#' @param reg_net An `igraph` object.
#' @param aracne A data frame.
#' @param MI_threshold Numeric.
#'
#' @return An `igraph` object.
#' @importFrom dplyr filter select mutate bind_rows
#' @importFrom igraph V as_data_frame graph_from_data_frame
#' @export
add_ARACNe <- function(reg_net, aracne, MI_threshold = NULL) {
  # Keep only ARACNe edges where both nodes exist in the network
  s_net <- dplyr::filter(aracne,
                         Regulator %in% igraph::V(reg_net)$name &
                         Target %in% igraph::V(reg_net)$name)

  # Apply MI threshold if provided
  if (!is.null(MI_threshold) && "MI" %in% colnames(s_net)) {
    s_net <- dplyr::filter(s_net, MI >= MI_threshold)
  }

  # Convert existing network to edge data frame
  edf <- igraph::as_data_frame(reg_net, what = "edges")

  # Add ARACNe edges
  edf <- dplyr::bind_rows(edf,
                          s_net |>
                            dplyr::select(from = Regulator, to = Target) |>
                            dplyr::mutate(source = "ARACNe"))

  # Build new igraph object
  nodes <- data.frame(name = igraph::V(reg_net)$name)
  new_net <- igraph::graph_from_data_frame(edf, vertices = nodes, directed = FALSE)

  return(new_net)
}

#' Add ARACNe edges to a regulatory network using MI threshold
#'
#' @param reg_net An igraph object.
#' @param dta A list.
#' @param MI_threshold Numeric or NULL.
#'
#' @return An igraph object.
#' @importFrom dplyr filter select mutate bind_rows rename
#' @importFrom igraph V as_data_frame graph_from_data_frame
#' @export
add_ARACNe_by_MI <- function(reg_net, dta, MI_threshold = NULL) {
  # Filter ARACNe edges to keep nodes present in nodes_all
  aracne_edges <- dta$aracne |>
    dplyr::filter(Regulator %in% dta$nodes_all$gene & Target %in% dta$nodes_all$gene)

  # Keep only edges where both nodes are in the current network
  aracne_edges <- aracne_edges |>
    dplyr::filter(Regulator %in% igraph::V(reg_net)$name & Target %in% igraph::V(reg_net)$name)

  # Apply MI threshold if specified
  if (!is.null(MI_threshold)) {
    aracne_edges <- aracne_edges |>
      dplyr::filter(MI >= MI_threshold)
  }

  # Convert reg_net to edge data frame
  reg_df <- igraph::as_data_frame(reg_net, what = "edges")

  # Add ARACNe edges
  reg_df <- dplyr::bind_rows(
    reg_df,
    aracne_edges |>
      dplyr::select(from = Regulator, to = Target) |>
      dplyr::mutate(source = "ARACNe")
  )

  # Build updated graph
  updated_net <- igraph::graph_from_data_frame(
    reg_df,
    vertices = dta$nodes_all |>
      dplyr::rename(name = gene) |>
      dplyr::filter(name %in% unique(c(reg_df$from, reg_df$to))),
    directed = FALSE
  )

  return(updated_net)
}

#' Save HotNet network files
#'
#' @param obj A list returned by `igraph_to_hotnet()`.
#' @param network_name Character string.
#' @param out_dir Character string.
#'
#' @return Invisibly returns a list of file paths.
#' @importFrom readr write_tsv
#' @export
save_hotnet_network_files <- function(obj, network_name, out_dir = "./data") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  index_file <- file.path(out_dir, paste0(network_name, "_index_gene.tsv"))
  edge_file  <- file.path(out_dir, paste0(network_name, "_edge_list.tsv"))
  edge_full_file <- file.path(out_dir, paste0(network_name, "_edge_list_full.tsv"))

  readr::write_tsv(obj$index, file = index_file, col_names = FALSE)
  readr::write_tsv(obj$edge,  file = edge_file,  col_names = FALSE)
  readr::write_tsv(obj$edge_full, file = edge_full_file, col_names = TRUE)

  message("HotNet files saved:\n", index_file, "\n", edge_file)

  invisible(list(index_file = index_file, edge_file = edge_file, edge_full_file = edge_full_file))
}

#' Select subgraphs containing specified nodes
#'
#' @param ig An `igraph` object.
#' @param nodes A character vector of node names.
#'
#' @return An `igraph` object.
#' @importFrom igraph V components induced_subgraph
#' @export
select_components_by_nodes <- function(ig, nodes) {
  # Ensure nodes exist in the graph
  nodes <- intersect(nodes, igraph::V(ig)$name)
  if (length(nodes) == 0) {
    stop("None of the specified nodes exist in the graph")
  }

  # Compute components
  comps <- igraph::components(ig)

  # Components that contain any of the specified nodes
  keep_components <- unique(comps$membership[match(nodes, igraph::V(ig)$name)])

  # Vertices to keep
  keep_vertices <- igraph::V(ig)[comps$membership %in% keep_components]

  # Return the induced subgraph
  igraph::induced_subgraph(ig, vids = keep_vertices)
}
