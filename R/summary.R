#' Load hotnet network input files
#'
#' @param net_prefix String specifying the network prefix (default: "network_1")
#' @param data_dir Path to data directory (required)
#'
#' @return A data frame containing network edges with columns from and to nodes
#' @importFrom readr read_tsv
#' @importFrom dplyr left_join rename select
#' @export
load_hotnet_input <- function(net_prefix = "network_1", data_dir) {
    if(is.null(data_dir)) {
        stop("data_dir parameter is required")
    }

    # Load index file
    idx <- readr::read_tsv(
        file.path(data_dir, paste0(net_prefix, "_index_gene.tsv")),
        col_names = c("idx", "node"),
        show_col_types = FALSE
    )

    # Load edge file
    edge <- readr::read_tsv(
        file.path(data_dir, paste0(net_prefix, "_edge_list.tsv")),
        col_names = c("id1", "id2"),
        show_col_types = FALSE
    )

    # Create and return network data frame as from-to format
    edge <- edge |>
        dplyr::left_join(idx, by = c("id1" = "idx")) |>
        dplyr::rename(from = "node") |>
        dplyr::left_join(idx, by = c("id2" = "idx")) |>
        dplyr::rename(to = "node") |>
        dplyr::select("from", "to")

    # load edge file with full info if exists
    edge_full_path <- file.path(data_dir, paste0(net_prefix, "_edge_list_full.tsv"))
    if(file.exists(edge_full_path)) {
        edge_full <- readr::read_tsv(
            edge_full_path,
            show_col_types = FALSE
        ) |> dplyr::select(-dplyr::any_of(c("id1", "id2")))
    }

    if(exists("edge_full")) {
        ### check consistency between edge and edge_full
        if(!all(edge$from == edge_full$from & edge$to == edge_full$to)) {
            stop("Inconsistency between edge and edge_full files")
        }
        return(edge_full)
    } else {
        return(edge)
    }
}

#' Load hotnet results from clusters file
#'
#' @param network String specifying the network name (default: "network_1")
#' @param result_dir Path to results directory (required)
#' @param exclude_single Logical indicating whether to exclude single-node clusters (default: FALSE)
#'
#' @return A list.
#' @importFrom stringr str_extract
#' @export
load_hotnet_results <- function(network = "network_1", result_dir, exclude_single = FALSE) {
    if(is.null(result_dir)) {
        stop("result_dir parameter is required")
    }

    # Find file starting with clusters_{network}
    files <- list.files(result_dir,
                       pattern = paste0("^clusters_", network, ".*\\.tsv$"),
                       full.names = TRUE)

    if(length(files) == 0) {
        stop("No matching results file found in: ", result_dir)
    }

    # Use first matching file if multiple exist
    filepath <- files[1]

    # Read all lines
    lines <- readLines(filepath)

    # Extract p-value from line containing "p-value:"
    pvalue_line <- lines[grep("p-value:", lines)[1]]
    pvalue <- as.numeric(stringr::str_extract(pvalue_line, "\\d+\\.\\d+$"))

    # Find first non-empty line after p-value line
    start_idx <- grep("p-value:", lines)[1] + 2

    # Read clusters starting from that line
    clusters <- lines[start_idx:length(lines)] |>
        # Split each line into vector of nodes
        strsplit(split = "\t") |>
        # Remove any empty elements
        lapply(function(x) x[x != ""])

    # Handle single-node clusters if requested
    if(exclude_single) {
        # Store single nodes
        singles <- clusters[lengths(clusters) == 1] |> unlist()
        # Keep only multi-node clusters
        clusters <- clusters[lengths(clusters) > 1]

        return(list(
            pvalue = pvalue,
            clusters = clusters,
            singles = singles
        ))
    }

    # Return results
    list(
        pvalue = pvalue,
        clusters = clusters
    )
}

#' Convert hotnet clusters list to a data frame
#'
#' @param res_net A list.
#'
#' @return A data frame.
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @export
hotnet_cluster_df <- function(res_net) {
    # Check input
    if (!("clusters" %in% names(res_net))) {
        stop("Input must contain a 'clusters' element")
    }

    # Create data frame from clusters list
    purrr::map_dfr(seq_along(res_net$clusters), function(i) {
        tibble::tibble(
            cluster = i,
            node = res_net$clusters[[i]]
        )
    })
}

#' Create igraph object from hotnet results with cluster-based node coloring and shapes
#'
#' @param edge_tbl Data frame of edges.
#' @param nodes_tbl Data frame of nodes.
#' @param hotnet Hotnet results.
#'
#' @return An igraph object.
#' @importFrom igraph graph_from_data_frame V
#' @importFrom dplyr filter left_join mutate coalesce if_else case_when
#' @importFrom stats setNames
#' @importFrom grDevices rainbow
#' @export
create_igraph_hotnet <- function(edge_tbl, nodes_tbl, hotnet) {
    # Get unique nodes from edge table
    nodes_in_edges <- unique(c(edge_tbl$from, edge_tbl$to))

    # Convert hotnet clusters to data frame
    clusters_df <- hotnet_cluster_df(hotnet)

    # Create color palette for clusters
    n_clusters <- if(nrow(clusters_df) > 0) max(clusters_df$cluster) else 0
    cluster_colors <- stats::setNames(
        grDevices::rainbow(n_clusters),
        1:n_clusters
    )

    # Add cluster information and shape to nodes table, filtering for nodes in edges
    nodes_with_clusters <- nodes_tbl |>
        dplyr::filter(gene %in% nodes_in_edges) |>
        dplyr::left_join(clusters_df, by = c("gene" = "node")) |>
        dplyr::mutate(
            cluster = dplyr::coalesce(cluster, 0),
            color = dplyr::if_else(
                cluster == 0,
                "gray",
                cluster_colors[as.character(cluster)]
            ),
            # Define shape based on source and node_type
            shape = dplyr::case_when(
                node_type == "Enhancer" & source == "JEME" ~ "square",
                node_type == "Enhancer" & source == "HiC" ~ "csquare",
                node_type == "Promoter" ~ "circle",
                TRUE ~ "circle"
            )
        )

    # Create igraph object
    g <- igraph::graph_from_data_frame(
        d = edge_tbl,
        vertices = nodes_with_clusters,
        directed = FALSE
    )

    # Set vertex attributes
    igraph::V(g)$color <- igraph::V(g)$color
    igraph::V(g)$shape <- igraph::V(g)$shape

    return(g)
}

#' Filter igraph object to keep only components with clustered nodes
#'
#' @param ig An igraph object.
#'
#' @return A filtered igraph object.
#' @importFrom igraph vertex_attr_names components vertex_attr induced_subgraph
#' @export
filter_igraph_hotnet <- function(ig) {
    # Check if cluster attribute exists
    if (!"cluster" %in% igraph::vertex_attr_names(ig)) {
        stop("Vertex attribute 'cluster' not found in the graph")
    }

    # Get components
    comps <- igraph::components(ig)

    # For each component, check if it contains any clustered nodes
    keep_components <- sapply(unique(comps$membership), function(cid) {
        # Get vertices in this component
        component_vertices <- which(comps$membership == cid)
        # Get their cluster values
        cluster_values <- igraph::vertex_attr(ig, "cluster")[component_vertices]
        # Return TRUE if any node has cluster != 0
        any(cluster_values != 0)
    })

    # Get vertices to keep
    vertices_to_keep <- which(comps$membership %in%
                            which(keep_components))

    # Create subgraph
    filtered_graph <- igraph::induced_subgraph(ig, vertices_to_keep)

    return(filtered_graph)
}

#' Plot each connected component of an igraph object to a multi-page PDF
#'
#' @param ig An igraph object.
#' @param score_attr Character string.
#' @param pdf_file Character string.
#' @param node_names Optional vector or logical.
#'
#' @return No return value.
#' @importFrom igraph vertex_attr_names components induced_subgraph vertex_attr V vcount layout_with_fr
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @export
plot_subgraphs_to_pdf2 <- function(ig, score_attr = "score", pdf_file = "subgraphs.pdf", node_names = FALSE) {
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
        # Node size proportional to score
        node_score <- igraph::vertex_attr(sg, score_attr)
        node_size <- 5 + 10 * (node_score / max(node_score, na.rm = TRUE))

        # Determine vertex labels
        vertex_labels <- if (isFALSE(node_names)) {
            NA
        } else if (isTRUE(node_names)) {
            igraph::V(sg)$name
        } else {
            ifelse(igraph::V(sg)$name %in% node_names,
                   igraph::V(sg)$name,
                   NA)
        }

        graphics::plot(
            sg,
            vertex.size = node_size,
            vertex.label = vertex_labels,
            vertex.shape = igraph::V(sg)$shape,
            vertex.color = igraph::V(sg)$color,
            vertex.label.cex = 0.6,
            vertex.label.color = "black",
            main = paste("Component with", igraph::vcount(sg), "nodes"),
            edge.color = "gray50",
            layout = igraph::layout_with_fr
        )
    }

    # Close PDF device
    grDevices::dev.off()

    message("Subgraphs saved to ", pdf_file)
}
