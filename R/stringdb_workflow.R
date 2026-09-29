stringdb_file_paths <- function(cache_dir, version, species) {
  prefix <- paste0(as.integer(species), ".protein.")
  suffix <- paste0(".v", as.character(version), ".txt.gz")
  c(
    aliases = file.path(cache_dir, paste0(prefix, "aliases", suffix)),
    info = file.path(cache_dir, paste0(prefix, "info", suffix)),
    links = file.path(cache_dir, paste0(prefix, "links", suffix))
  )
}

#' Download and cache a versioned STRING protein network
#'
#' This is an explicit data-provisioning step. Analysis execution requires all
#' three versioned cache files to exist and never downloads them implicitly.
#'
#' @param cache_dir Destination directory. Defaults to `.cache/STRINGdb` below
#'   the current project cache.
#' @param version STRING database release.
#' @param species NCBI taxonomy identifier.
#' @param score_threshold Combined-score threshold used to initialize STRINGdb.
#' @param overwrite Remove and reacquire the three versioned cache files.
#' @return A data frame containing each required path and its status.
#' @export
download_stringdb_data <- function(cache_dir = NULL, version = "12.0",
                                   species = 9606L,
                                   score_threshold = 700,
                                   overwrite = FALSE) {
  assert_flag(overwrite, "overwrite")
  assert_positive_integer(species, "species")
  assert_scalar(version, "version", "character")
  assert_scalar(score_threshold, "score_threshold")
  if (!is.numeric(score_threshold) || !is.finite(score_threshold) ||
      score_threshold < 0 || score_threshold > 1000) {
    stop("score_threshold must be between 0 and 1000", call. = FALSE)
  }
  if (is.null(cache_dir)) {
    cache_dir <- file.path(get_cache_dir(), "STRINGdb")
  }
  cache_dir <- normalizePath(cache_dir, mustWork = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- stringdb_file_paths(cache_dir, version, species)
  if (isTRUE(overwrite)) unlink(paths[file.exists(paths)])
  before <- file.exists(paths)
  if (!all(before)) {
    STRINGdb::STRINGdb$new(
      version = version, species = as.integer(species),
      score_threshold = score_threshold, input_directory = cache_dir
    )
  }
  available <- file.exists(paths) & file.info(paths)$size > 0
  if (!all(available)) {
    stop(
      "STRINGdb setup did not produce: ",
      paste(paths[!available], collapse = ", "), call. = FALSE
    )
  }
  data.frame(
    artifact = names(paths), path = unname(paths),
    status = ifelse(before & !overwrite, "existing", "downloaded"),
    stringsAsFactors = FALSE
  )
}

new_analysis_stringdb <- function(config) {
  settings <- stringdb_settings(config)
  paths <- stringdb_cache_paths(config)
  missing <- paths[!file.exists(paths) | file.info(paths)$size == 0]
  if (length(missing)) {
    stop(
      "Missing STRINGdb cache file(s): ", paste(missing, collapse = ", "),
      ". Run download_stringdb_data() explicitly before analysis.",
      call. = FALSE
    )
  }
  new_cached_stringdb_client(paths, settings$score_threshold)
}

stringdb_mapping_key <- function(x) {
  toupper(iconv(as.character(x), from = "WINDOWS-1252", to = "UTF-8"))
}

new_cached_stringdb_client <- function(paths, score_threshold) {
  info <- readr::read_tsv(
    paths[["stringdb_info"]], show_col_types = FALSE, progress = FALSE,
    col_select = 1:2
  )
  aliases <- readr::read_tsv(
    paths[["stringdb_aliases"]], show_col_types = FALSE, progress = FALSE,
    col_select = 1:2
  )
  names(info) <- c("STRING_id", "preferred_name")
  names(aliases) <- c("STRING_id", "alias")
  preferred <- dplyr::bind_rows(
    dplyr::transmute(info, STRING_id = as.character(.data$STRING_id),
                     alias = as.character(.data$preferred_name)),
    dplyr::transmute(info, STRING_id = as.character(.data$STRING_id),
                     alias = as.character(.data$STRING_id))
  )
  preferred_keys <- unique(stringdb_mapping_key(preferred$alias))
  lookup <- aliases |>
    dplyr::transmute(
      STRING_id = as.character(.data$STRING_id),
      alias = as.character(.data$alias),
      mapping_key = stringdb_mapping_key(.data$alias)
    ) |>
    dplyr::filter(!.data$mapping_key %in% preferred_keys) |>
    dplyr::bind_rows(
      dplyr::mutate(preferred,
                    mapping_key = stringdb_mapping_key(.data$alias))
    ) |>
    dplyr::select("mapping_key", "STRING_id") |>
    dplyr::distinct()
  links <- NULL

  map <- function(my_data_frame, my_data_frame_id_col_names,
                  takeFirst = TRUE, removeUnmappedRows = FALSE,
                  quiet = FALSE) {
    if (length(my_data_frame_id_col_names) != 1L) {
      stop("Cached STRINGdb mapping currently accepts one identifier column",
           call. = FALSE)
    }
    id_column <- my_data_frame_id_col_names[[1L]]
    if (!id_column %in% names(my_data_frame)) {
      stop("Missing STRINGdb mapping column: ", id_column, call. = FALSE)
    }
    input <- dplyr::mutate(
      my_data_frame, .stringdb_row = dplyr::row_number(),
      .stringdb_key = stringdb_mapping_key(.data[[id_column]])
    )
    mapped <- dplyr::left_join(
      input, lookup, by = c(".stringdb_key" = "mapping_key"),
      relationship = "many-to-many"
    ) |>
      dplyr::arrange(.data$.stringdb_row) |>
      dplyr::select(-".stringdb_row", -".stringdb_key")
    if (isTRUE(removeUnmappedRows)) {
      mapped <- dplyr::filter(mapped, !is.na(.data$STRING_id))
    }
    mapped
  }

  get_interactions <- function(string_ids) {
    if (is.null(links)) {
      links <<- readr::read_table(
        paths[["stringdb_links"]], show_col_types = FALSE, progress = FALSE,
        col_types = readr::cols(
          protein1 = readr::col_character(),
          protein2 = readr::col_character(),
          combined_score = readr::col_double()
        )
      ) |>
        dplyr::filter(.data$combined_score >= score_threshold)
    }
    links |>
      dplyr::filter(
        .data$protein1 %in% string_ids,
        .data$protein2 %in% string_ids
      ) |>
      dplyr::transmute(
        from = .data$protein1, to = .data$protein2,
        combined_score = .data$combined_score
      )
  }

  list(map = map, get_interactions = get_interactions)
}

empty_stringdb_result <- function(graph, eligible_genes = character()) {
  list(
    graph = graph,
    mapping = data.frame(
      gene = eligible_genes, STRING_id = NA_character_,
      mapping_status = "unmapped", selected_for_string_id = FALSE,
      stringsAsFactors = FALSE
    ),
    conflicts = data.frame(
      STRING_id = character(), gene_count = integer(), genes = character(),
      selected_gene = character(), stringsAsFactors = FALSE
    ),
    edges = data.frame(
      Regulator = character(), Target = character(), source = character(),
      combined_score = numeric(), Regulator_STRING_id = character(),
      Target_STRING_id = character(), overlaps_base = logical(),
      stringsAsFactors = FALSE
    ),
    summary = data.frame(
      eligible_gene_nodes = length(eligible_genes), mapped_gene_rows = 0L,
      mapped_string_ids = 0L, mapping_conflicts = 0L,
      interactions_at_threshold = 0L, self_loops_removed = 0L,
      duplicate_gene_pairs_removed = 0L, overlaps_existing_base = 0L,
      string_edges_added = 0L, stringsAsFactors = FALSE
    )
  )
}

prepare_stringdb_augmentation <- function(graph, stringdb_obj,
                                          score_threshold = 700,
                                          mapping_policy =
                                            "legacy_shortest_gene_per_string_id") {
  if (!identical(mapping_policy, "legacy_shortest_gene_per_string_id")) {
    stop("Unsupported STRINGdb mapping policy: ", mapping_policy, call. = FALSE)
  }
  vertices <- igraph::as_data_frame(graph, what = "vertices")
  if (!all(c("name", "node_type") %in% names(vertices))) {
    stop("STRINGdb augmentation requires name and node_type vertex attributes",
         call. = FALSE)
  }
  genes <- sort(unique(as.character(
    vertices$name[vertices$node_type == "promoter_gene"]
  )))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) return(empty_stringdb_result(graph))

  mapped_raw <- stringdb_obj$map(
    data.frame(gene = genes, stringsAsFactors = FALSE),
    "gene", removeUnmappedRows = FALSE
  )
  if (!all(c("gene", "STRING_id") %in% names(mapped_raw))) {
    stop("STRINGdb mapping did not return gene and STRING_id columns",
         call. = FALSE)
  }
  mapping <- mapped_raw |>
    dplyr::transmute(
      gene = as.character(.data$gene),
      STRING_id = as.character(.data$STRING_id),
      mapping_status = dplyr::if_else(
        is.na(.data$STRING_id) | !nzchar(.data$STRING_id),
        "unmapped", "mapped"
      )
    ) |>
    dplyr::distinct()
  mapped <- mapping |>
    dplyr::filter(.data$mapping_status == "mapped") |>
    dplyr::arrange(.data$STRING_id, nchar(.data$gene), .data$gene) |>
    dplyr::mutate(
      selected_for_string_id = !duplicated(.data$STRING_id)
    )
  mapping <- mapping |>
    dplyr::left_join(
      dplyr::select(mapped, "gene", "STRING_id", "selected_for_string_id"),
      by = c("gene", "STRING_id")
    ) |>
    dplyr::mutate(
      selected_for_string_id = dplyr::coalesce(
        .data$selected_for_string_id, FALSE
      )
    ) |>
    dplyr::arrange(.data$gene, .data$STRING_id)
  selected <- dplyr::filter(mapped, .data$selected_for_string_id)
  conflicts <- selected |>
    dplyr::select("STRING_id", selected_gene = "gene") |>
    dplyr::right_join(
      mapped |>
        dplyr::summarise(
          gene_count = dplyr::n_distinct(.data$gene),
          genes = paste(sort(unique(.data$gene)), collapse = ";"),
          .by = "STRING_id"
        ) |>
        dplyr::filter(.data$gene_count > 1L),
      by = "STRING_id"
    ) |>
    dplyr::select("STRING_id", "gene_count", "genes", "selected_gene") |>
    dplyr::arrange(.data$STRING_id)
  if (!nrow(selected)) {
    result <- empty_stringdb_result(graph, genes)
    result$mapping <- mapping
    result$conflicts <- conflicts
    return(result)
  }

  interactions <- stringdb_obj$get_interactions(unique(selected$STRING_id))
  required <- c("from", "to", "combined_score")
  missing <- setdiff(required, names(interactions))
  if (length(missing)) {
    stop(
      "STRINGdb interactions are missing: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  interactions <- interactions |>
    dplyr::transmute(
      from_string = as.character(.data$from),
      to_string = as.character(.data$to),
      combined_score = as.numeric(.data$combined_score)
    ) |>
    dplyr::filter(
      is.finite(.data$combined_score),
      .data$combined_score >= score_threshold
    ) |>
    dplyr::distinct()
  threshold_rows <- nrow(interactions)
  gene_edges <- interactions |>
    dplyr::inner_join(
      dplyr::select(selected, from_string = "STRING_id", from_gene = "gene"),
      by = "from_string", relationship = "many-to-one"
    ) |>
    dplyr::inner_join(
      dplyr::select(selected, to_string = "STRING_id", to_gene = "gene"),
      by = "to_string", relationship = "many-to-one"
    )
  self_loops <- sum(gene_edges$from_gene == gene_edges$to_gene)
  gene_edges <- gene_edges |>
    dplyr::filter(.data$from_gene != .data$to_gene) |>
    dplyr::mutate(
      forward = .data$from_gene < .data$to_gene,
      Regulator = ifelse(.data$forward, .data$from_gene, .data$to_gene),
      Target = ifelse(.data$forward, .data$to_gene, .data$from_gene),
      Regulator_STRING_id = ifelse(
        .data$forward, .data$from_string, .data$to_string
      ),
      Target_STRING_id = ifelse(
        .data$forward, .data$to_string, .data$from_string
      )
    ) |>
    dplyr::arrange(
      .data$Regulator, .data$Target, dplyr::desc(.data$combined_score),
      .data$Regulator_STRING_id, .data$Target_STRING_id
    )
  before_gene_dedup <- nrow(gene_edges)
  gene_edges <- gene_edges |>
    dplyr::distinct(.data$Regulator, .data$Target, .keep_all = TRUE)
  duplicate_pairs <- before_gene_dedup - nrow(gene_edges)

  base <- igraph::as_data_frame(graph, what = "edges")
  base_from <- as.character(base$from)
  base_to <- as.character(base$to)
  base_keys <- paste(pmin(base_from, base_to), pmax(base_from, base_to), sep = "\r")
  string_keys <- paste(gene_edges$Regulator, gene_edges$Target, sep = "\r")
  gene_edges <- gene_edges |>
    dplyr::mutate(
      source = "STRINGdb",
      overlaps_base = string_keys %in% base_keys
    ) |>
    dplyr::select(
      "Regulator", "Target", "source", "combined_score",
      "Regulator_STRING_id", "Target_STRING_id", "overlaps_base"
    )
  add <- dplyr::filter(gene_edges, !.data$overlaps_base) |>
    dplyr::transmute(
      from = .data$Regulator, to = .data$Target, source = .data$source,
      weight = .data$combined_score,
      string_score = .data$combined_score,
      from_string = .data$Regulator_STRING_id,
      to_string = .data$Target_STRING_id
    )
  combined <- dplyr::bind_rows(base, add)
  augmented <- igraph::graph_from_data_frame(
    combined, vertices = vertices, directed = FALSE
  )
  list(
    graph = augmented,
    mapping = mapping,
    conflicts = conflicts,
    edges = gene_edges,
    summary = data.frame(
      eligible_gene_nodes = length(genes),
      mapped_gene_rows = nrow(dplyr::filter(mapping, .data$mapping_status == "mapped")),
      mapped_string_ids = dplyr::n_distinct(selected$STRING_id),
      mapping_conflicts = nrow(conflicts),
      interactions_at_threshold = threshold_rows,
      self_loops_removed = self_loops,
      duplicate_gene_pairs_removed = duplicate_pairs,
      overlaps_existing_base = sum(gene_edges$overlaps_base),
      string_edges_added = nrow(add), stringsAsFactors = FALSE
    )
  )
}
