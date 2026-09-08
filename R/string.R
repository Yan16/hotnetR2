#' Add STRINGdb Edges to an igraph Object
#'
#' This function takes an igraph object and a STRINGdb object, retrieves the edges from STRINGdb
#' that correspond to the nodes in the igraph object, and combines them with the existing edges
#' in the igraph object.
#'
#' @param ig An igraph object representing the existing network.
#' @param stringdb_obj A STRINGdb object.
#' @param simplify Logical, if TRUE, remove duplicated edges and loops (default = FALSE).
#'
#' @return An igraph object.
#' @importFrom igraph as_data_frame graph_from_data_frame simplify
#' @importFrom dplyr select rename right_join filter group_by n_distinct slice_min mutate relocate bind_rows any_of
#' @importFrom stats na.omit
#' @export
add_STRINGdb_edges <- function(ig, stringdb_obj, simplify=FALSE) {

    ##### 1. data from igraph object
    ##   get the edge table
    edge_list <- igraph::as_data_frame(ig, what = "edges")

    ## get the node table
    node_list <- igraph::as_data_frame(ig, what = "vertices") |>
      dplyr::select(-dplyr::any_of(c("cluster", "color", "shape")))

    #### 2. add STRINGdb data to node list
    ## Add STRINGdb mapping to node_list
    node_update <- stringdb_obj$map(as.data.frame(node_list), "name", removeUnmappedRows = FALSE)

    ## only those mapped to STING ids
    node_maped <- node_update |> dplyr::select(name, STRING_id) |> dplyr::filter(!is.na(STRING_id))

    # Find IDs that map to multiple names
    conflicts <- node_maped |>
      dplyr::group_by(STRING_id) |>
      dplyr::filter(dplyr::n_distinct(name) > 1) |>
      dplyr::ungroup()

    # Display the conflicting rows
    if(nrow(conflicts) > 0) {
        warning("Multiple node/gene mapped to single STRING_id")
        print(conflicts)
        print("The row with shortest gene name will be kept")

        node_maped <- node_maped |>
          dplyr::group_by(STRING_id) |>
          # keep the row with the shortest name string
          dplyr::slice_min(nchar(name), n = 1, with_ties = FALSE) |>
          dplyr::ungroup()
    }

    maped_id <- stats::na.omit(node_maped$STRING_id)


    #### 3. identify edges in STRINGdb edges that link nodes in igraph object
    edges <- stringdb_obj$get_interactions(maped_id) |> dplyr::distinct()
    # Prepare edge data frame for igraph, by mapping STRING IDs back to gene names
    # and adding source column
    edge_df <- edges |>
        dplyr::rename(from_string = from, to_string = to, string_score = combined_score) |>
        dplyr::right_join(
          node_maped |> dplyr::select(STRING_id, from = name),
          by = c("from_string" = "STRING_id")
        ) |>
        dplyr::filter(!is.na(from)) |>
        dplyr::right_join(
          node_maped |> dplyr::select(STRING_id, to = name),
          by = c("to_string" = "STRING_id")
        ) |>
        dplyr::mutate(source = "STRINGdb") |>
        dplyr::relocate(from, to, source, dplyr::everything()) |>
        dplyr::filter(!is.na(from) & !is.na(to))

  # Combine with existing network edges
  combined_edges <- igraph::as_data_frame(ig, "edges") |>
    dplyr::bind_rows(edge_df)

  # Build vertex table from
  vertices <- node_update

  # Create new igraph object
  out <- igraph::graph_from_data_frame(combined_edges, vertices = vertices, directed = FALSE)
  if(simplify) {
    out <- igraph::simplify(out, remove.multiple = TRUE, remove.loops = TRUE)
  }
  out
}
