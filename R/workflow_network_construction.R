# Construct the four standard h-HotNet input networks inside analysis1.

require_network_packages <- function() {
  packages <- c("dplyr", "readr", "stringr", "tidyr", "igraph")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Network construction requires: ", paste(missing, collapse = ", "), call. = FALSE)
}

network_data_file <- function(config, name) {
  path <- file.path(analysis_paths(config)[["hhotnet_data"]], name)
  root <- analysis_paths(config)[["root"]]
  if (!startsWith(normalizePath(path, mustWork = FALSE), paste0(root, .Platform$file.sep))) {
    stop("Refusing network output outside analysis directory", call. = FALSE)
  }
  path
}

network_summary_file <- function(config, name) {
  path <- file.path(analysis_paths(config)[["hhotnet_summary"]], name)
  root <- analysis_paths(config)[["root"]]
  if (!startsWith(normalizePath(path, mustWork = FALSE), paste0(root, .Platform$file.sep))) {
    stop("Refusing network output outside analysis directory", call. = FALSE)
  }
  path
}

normalise_network_node_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^chr", "chr", x, ignore.case = TRUE)
  gsub("[.:-]", "_", x)
}

score_from_pvalue <- function(p, cap = NULL) {
  p <- suppressWarnings(as.numeric(p))
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("LDAK p-values must lie in [0, 1]", call. = FALSE)
  score <- -log10(p)
  finite <- score[is.finite(score)]
  if (is.null(cap)) cap <- if (length(finite)) max(finite) + 1 else 1
  if (!is.numeric(cap) || length(cap) != 1L || !is.finite(cap) || cap <= 0) {
    stop("scores.infinite_score_cap must be null or a positive finite number", call. = FALSE)
  }
  score[is.infinite(score) & score > 0] <- cap
  score
}

extract_gene_nodes <- function(jeme_nodes, hic_nodes, aracne) {
  jeme_promoters <- jeme_nodes$promoters
  hic_promoters <- hic_nodes$promoters
  jeme_gene <- if (!is.null(jeme_promoters) && nrow(jeme_promoters)) jeme_promoters$promoter %||% jeme_promoters$Promoter else character()
  hic_gene <- if (!is.null(hic_promoters) && nrow(hic_promoters)) hic_promoters$Promoter %||% hic_promoters$promoter else character()
  aracne_gene <- if (!is.null(aracne) && nrow(aracne)) unique(c(aracne$Regulator, aracne$Target)) else character()
  dplyr::bind_rows(
    data.frame(gene = jeme_gene, source = rep("JEME_promoter", length(jeme_gene))),
    data.frame(gene = hic_gene, source = rep("HiC_promoter", length(hic_gene))),
    data.frame(gene = aracne_gene, source = rep("ARACNe", length(aracne_gene)))
  ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene)) |>
    dplyr::summarise(sources = paste(sort(unique(.data$source)), collapse = ";"), .by = "gene") |>
    dplyr::arrange(.data$gene)
}

refresh_network_nodes <- function(data) {
  data$nodes$jeme <- extract_jeme_nodes(data$jeme)
  data$nodes$hic <- extract_hic_nodes(data$hic)
  data$nodes$all <- extract_gene_nodes(data$nodes$jeme, data$nodes$hic, data$aracne)
  data
}

load_network_sources <- function(config) {
  jeme <- config$regulatory$jeme
  if (is.null(jeme$key_column) || is.null(jeme$value)) {
    stop("regulatory.jeme.key_column/value are required", call. = FALSE)
  }
  jeme_data <- get_jeme(
    key_column = jeme$key_column, value = jeme$value,
    method = jeme$method %||% "lasso", simplified = jeme$simplified %||% TRUE,
    cache_dir = config$cache_dir
  )
  hic_data <- load_analysis_hic(config)
  if (isTRUE(config$aracne$enabled %||% TRUE)) {
    aracne_file <- resolve_config_path(config$aracne$file, config$project_root)
    if (!file.exists(aracne_file)) stop("Missing ARACNe network: ", aracne_file, call. = FALSE)
    aracne_data <- read_ARACNe(
      infile = aracne_file, MI_threshold = config$aracne$MI_threshold %||% 0,
      pvalue_threshold = config$aracne$pvalue_threshold %||% 1e-6,
      keep_genes = config$aracne$keep_genes %||% NULL
    )
    aracne_harmonization <- harmonize_aracne_endpoints(
      aracne_data,
      config$aracne$harmonization %||% NULL,
      config$project_root
    )
  } else {
    aracne_harmonization <- list(
      data = data.frame(
        Regulator = character(), Target = character(), MI = numeric(),
        stringsAsFactors = FALSE
      ),
      mapping = data.frame(
        endpoint = character(), original = character(), canonical = character(),
        mapping_status = character(), mapping_source = character(),
        mapping_candidates = character(), stringsAsFactors = FALSE
      ),
      endpoint_summary = data.frame(
        endpoint = character(), mapping_status = character(),
        identifiers = integer(), occurrences = integer(), stringsAsFactors = FALSE
      ),
      edge_summary = data.frame(
        harmonization_enabled = FALSE,
        input_edges = 0L,
        harmonized_edges = 0L,
        dropped_edges = 0L,
        reason = "ARACNe disabled in analysis.yaml"
      )
    )
  }
  data <- list(
    jeme = jeme_data, hic = hic_data,
    aracne = aracne_harmonization$data,
    aracne_harmonization = aracne_harmonization,
    nodes = list()
  )
  refresh_network_nodes(data)
}

deduplicate_network_jeme <- function(data, enabled = FALSE) {
  before <- nrow(dplyr::bind_rows(data$jeme))
  if (isTRUE(enabled)) {
    data$jeme <- list(dplyr::bind_rows(data$jeme) |>
      dplyr::distinct(.data$enhancer, .data$promoter, .keep_all = TRUE))
    data <- refresh_network_nodes(data)
  }
  list(data = data, summary = data.frame(
    enabled = isTRUE(enabled), rows_before = before, rows_after = nrow(dplyr::bind_rows(data$jeme)),
    rows_removed = before - nrow(dplyr::bind_rows(data$jeme))
  ))
}

split_network_multi_promoters <- function(data, enabled = FALSE) {
  po_before <- nrow(data$hic$PO); pp_before <- nrow(data$hic$PP)
  if (isTRUE(enabled)) {
    data$hic$PO <- data$hic$PO |>
      dplyr::distinct(.data$Interacting_fragment, .data$Promoter, .keep_all = TRUE) |>
      tidyr::separate_rows("Promoter", sep = ";") |>
      dplyr::mutate(Promoter = stringr::str_trim(.data$Promoter))
    data$hic$PP <- data$hic$PP |>
      tidyr::separate_rows("Promoter1", sep = ";") |>
      dplyr::mutate(Promoter1 = stringr::str_trim(.data$Promoter1)) |>
      tidyr::separate_rows("Promoter2", sep = ";") |>
      dplyr::mutate(Promoter2 = stringr::str_trim(.data$Promoter2))
    data <- refresh_network_nodes(data)
  }
  list(data = data, summary = data.frame(
    enabled = isTRUE(enabled), hic_po_rows_before = po_before, hic_po_rows_after = nrow(data$hic$PO),
    hic_pp_rows_before = pp_before, hic_pp_rows_after = nrow(data$hic$PP)
  ))
}

harmonize_network_promoters <- function(data, settings, analysis_config) {
  if (!isTRUE(settings$enabled %||% FALSE)) {
    return(list(data = data, summary = data.frame(enabled = FALSE, jeme_changed = 0L, hic_po_changed = 0L, hic_pp_changed = 0L)))
  }
  changed_jeme <- changed_hic_po <- changed_hic_pp <- 0L
  jeme_mapping <- data.frame()
  if (isTRUE(settings$jeme %||% TRUE)) {
    mapping_file <- annotation_output_file(
      analysis_config, "jeme_ensg_harmonization.tsv.gz"
    )
    if (!file.exists(mapping_file)) {
      stop(
        "Missing all-tissue JEME ENSG mapping: ", mapping_file,
        ". Regenerate LDAK annotations before constructing networks.",
        call. = FALSE
      )
    }
    jeme_mapping <- readr::read_tsv(
      mapping_file, show_col_types = FALSE, progress = FALSE
    )
    input <- dplyr::bind_rows(data$jeme)
    output <- apply_all_tissue_jeme_ensg_map(input, jeme_mapping)
    changed_jeme <- sum(output$promoter != output$promoter_original, na.rm = TRUE)
    data$jeme <- list(output)
  }
  if (analysis_hic_enabled(analysis_config) && isTRUE(settings$hic %||% TRUE)) {
    alias_file <- settings$hic_alias_file %||% NULL
    if (is.null(alias_file)) stop("regulatory.promoter_harmonization.hic_alias_file is required when Hi-C harmonization is enabled", call. = FALSE)
    alias_file <- resolve_config_path(alias_file, analysis_config$project_root)
    if (!file.exists(alias_file)) stop("Missing Hi-C promoter alias mapping: ", alias_file, call. = FALSE)
    hic_result <- harmonize_hic_gene_symbols(data$hic, alias_file)
    data$hic <- hic_result$data
    changed_hic_po <- hic_result$summary$hic_po_changed
    changed_hic_pp <- hic_result$summary$hic_pp_changed
  }
  data <- refresh_network_nodes(data)
  list(
    data = data,
    jeme_mapping = jeme_mapping,
    summary = data.frame(
      enabled = TRUE,
      jeme_method = "all_tissues_ENSG_map_then_selected_tissue",
      jeme_map_ensg = nrow(jeme_mapping),
      jeme_changed = changed_jeme,
      hic_method = "gene_symbol_alias_mapping",
      hic_po_changed = changed_hic_po,
      hic_pp_changed = changed_hic_pp
    )
  )
}

classify_network_enhancer_edges <- function(data, config) {
  settings <- config$enhancer_classification %||% list(enabled = FALSE)
  if (!isTRUE(settings$enabled %||% FALSE)) {
    return(list(
      data = data,
      audit = data.frame(),
      summary = data.frame(enabled = FALSE)
    ))
  }
  class1 <- as.character(settings$class1_label %||% "enh_class1")
  class2 <- as.character(settings$class2_label %||% "enh_class2")
  class3 <- as.character(settings$class3_label %||% "enh_class3")
  drop_sources <- enhancer_drop_sources(config)
  classification_file <- annotation_output_file(config, "enhancer_classification.tsv.gz")
  if (!file.exists(classification_file)) {
    stop("Missing enhancer classification table: ", classification_file, call. = FALSE)
  }
  classes <- readr::read_tsv(
    classification_file, show_col_types = FALSE, progress = FALSE,
    col_select = c("name", "enhancer_class")
  )
  invalid_classes <- setdiff(unique(classes$enhancer_class), c(class1, class2, class3))
  if (length(invalid_classes) || anyDuplicated(classes$name)) {
    stop("Enhancer classification contains invalid classes or duplicate names", call. = FALSE)
  }

  jeme <- dplyr::bind_rows(data$jeme) |>
    dplyr::left_join(classes, by = c("enhancer" = "name"))
  hic_po <- data$hic$PO |>
    dplyr::left_join(classes, by = c("Interacting_fragment" = "name"))
  missing_jeme <- sum(is.na(jeme$enhancer_class))
  missing_hic <- sum(is.na(hic_po$enhancer_class))
  if (missing_jeme || missing_hic) {
    stop(
      "Enhancer class is missing for ", missing_jeme,
      " JEME and ", missing_hic, " HiC-PO edge(s)", call. = FALSE
    )
  }
  audit <- dplyr::bind_rows(
    dplyr::transmute(
      jeme, Regulator = .data$enhancer, Target = .data$promoter,
      source = "JEME", enhancer_class = .data$enhancer_class,
      edge_action = dplyr::if_else(.data$enhancer_class == class1 | !"JEME" %in% drop_sources, "KEEP", "DROP")
    ),
    dplyr::transmute(
      hic_po, Regulator = .data$Interacting_fragment, Target = .data$Promoter,
      source = "HiC_PO", enhancer_class = .data$enhancer_class,
      edge_action = dplyr::if_else(.data$enhancer_class == class1 | !"HiC_PO" %in% drop_sources, "KEEP", "DROP")
    )
  )
  data$jeme <- list(dplyr::filter(jeme, .data$enhancer_class == class1 | !"JEME" %in% drop_sources) |>
                      dplyr::select(-dplyr::all_of("enhancer_class")))
  data$hic$PO <- dplyr::filter(hic_po, .data$enhancer_class == class1 | !"HiC_PO" %in% drop_sources) |>
    dplyr::select(-dplyr::all_of("enhancer_class"))
  data <- refresh_network_nodes(data)
  summary <- dplyr::summarise(
    dplyr::group_by(audit, .data$source, .data$enhancer_class, .data$edge_action),
    edge_rows = dplyr::n(), unique_enhancers = dplyr::n_distinct(.data$Regulator),
    .groups = "drop"
  ) |>
    dplyr::mutate(enabled = TRUE, .before = 1L)
  list(data = data, audit = audit, summary = summary)
}

filter_analysis_network_ldak <- function(data, config) {
  paths <- ldak_summary_specs(config)
  enhancer_file <- paths$enhancer$output_file; promoter_file <- paths$promoter$output_file
  if (!file.exists(enhancer_file) || !file.exists(promoter_file)) stop("Standardized LDAK summaries are required before network construction", call. = FALSE)
  enhancer_nodes <- dplyr::bind_rows(
    dplyr::transmute(data$nodes$jeme$enhancers, node = .data$enhancer, node_sources = "JEME"),
    dplyr::transmute(data$nodes$hic$enhancers, node = .data$Interacting_fragment, node_sources = "HiC")
  ) |>
    dplyr::filter(!is.na(.data$node), nzchar(.data$node)) |>
    dplyr::mutate(node_key = normalise_network_node_id(.data$node)) |>
    dplyr::summarise(node = paste(sort(unique(.data$node)), collapse = ";"), node_sources = paste(sort(unique(.data$node_sources)), collapse = ";"), .by = "node_key")
  promoter_nodes <- data$nodes$all |>
    dplyr::transmute(node = .data$gene, node_sources = .data$sources, node_key = normalise_network_node_id(.data$gene)) |>
    dplyr::filter(!is.na(.data$node), nzchar(.data$node)) |>
    dplyr::distinct()
  enhancers <- readr::read_tsv(enhancer_file, show_col_types = FALSE, progress = FALSE) |>
    dplyr::mutate(node_key = normalise_network_node_id(.data$gene)) |>
    dplyr::inner_join(enhancer_nodes, by = "node_key") |>
    dplyr::select(-"node_key")
  promoters <- readr::read_tsv(promoter_file, show_col_types = FALSE, progress = FALSE) |>
    dplyr::mutate(node_key = normalise_network_node_id(.data$gene)) |>
    dplyr::inner_join(promoter_nodes, by = "node_key") |>
    dplyr::select(-"node_key")
  list(
    enhancers = enhancers, promoters = promoters,
    unmatched_enhancers = dplyr::anti_join(enhancer_nodes, dplyr::mutate(enhancers, node_key = normalise_network_node_id(.data$gene)), by = "node_key"),
    unmatched_promoters = dplyr::anti_join(promoter_nodes, dplyr::mutate(promoters, node_key = normalise_network_node_id(.data$gene)), by = "node_key")
  )
}

exclude_network_mhc <- function(nodes, config) {
  if (!isTRUE(config$enabled %||% FALSE)) {
    return(list(
      data = nodes,
      excluded = nodes[0, , drop = FALSE],
      region = mhc_exclusion_region(config)
    ))
  }
  required <- c("gene", "CHR", "START", "END")
  missing <- setdiff(required, names(nodes))
  if (length(missing)) {
    stop(
      "MHC exclusion requires columns: ",
      paste(required, collapse = ", "),
      "; missing: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  region <- mhc_exclusion_region(config)
  chromosome <- toupper(trimws(as.character(nodes$CHR)))
  chromosome <- sub("^CHR", "", chromosome)
  start <- suppressWarnings(as.numeric(nodes$START))
  end <- suppressWarnings(as.numeric(nodes$END))
  invalid <- is.na(chromosome) | !nzchar(chromosome) |
    !is.finite(start) | !is.finite(end) | start < 0 | start >= end
  if (any(invalid)) {
    stop(
      "MHC exclusion cannot guarantee complete removal because ",
      sum(invalid), " node interval(s) have missing or invalid coordinates",
      call. = FALSE
    )
  }
  region_start <- region$start - 1L
  region_end <- region$end
  overlaps <- chromosome == as.character(region$chr) &
    end > region_start & start < region_end
  excluded_ids <- unique(as.character(nodes$gene[overlaps]))
  excluded <- as.character(nodes$gene) %in% excluded_ids
  list(
    data = nodes[!excluded, , drop = FALSE],
    excluded = nodes[excluded, , drop = FALSE],
    region = region
  )
}

apply_network_self_loop_policy <- function(edges, preserve_self_loops = TRUE) {
  if (isTRUE(preserve_self_loops)) {
    edges
  } else {
    dplyr::filter(edges, .data$Regulator != .data$Target)
  }
}

filter_network_regulatory_edges <- function(edges, nodes, mode,
                                            preserve_self_loops = TRUE) {
  if (!mode %in% c("current", "legacy_filter_network_ldak")) stop("Unknown networks.regulatory_edge_filter: ", mode, call. = FALSE)
  candidates <- if (identical(mode, "current")) {
    edges |>
      dplyr::filter(!is.na(.data$Regulator), !is.na(.data$Target), nzchar(.data$Regulator), nzchar(.data$Target)) |>
      dplyr::distinct()
  } else edges
  candidates <- apply_network_self_loop_policy(candidates, preserve_self_loops)
  keep <- candidates$Regulator %in% nodes$gene & candidates$Target %in% nodes$gene
  list(kept = candidates[keep, , drop = FALSE], dropped = candidates[!keep, , drop = FALSE])
}

make_network_graph <- function(edges, nodes) {
  vertices <- nodes |>
    dplyr::filter(.data$gene %in% unique(c(edges$Regulator, edges$Target))) |>
    dplyr::rename(name = "gene") |>
    dplyr::relocate("name")
  missing <- setdiff(unique(c(edges$Regulator, edges$Target)), vertices$name)
  if (length(missing)) stop("Network edges refer to missing node(s): ", paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
  igraph::graph_from_data_frame(
    edges |> dplyr::rename(from = "Regulator", to = "Target"),
    vertices = vertices, directed = FALSE
  )
}

subset_network_components <- function(graph, threshold) {
  components <- igraph::components(graph)
  scores <- igraph::vertex_attr(graph, "score")
  keep_component <- tapply(scores, components$membership, max, na.rm = TRUE) >= threshold
  igraph::induced_subgraph(graph, vids = igraph::V(graph)[keep_component[as.character(components$membership)]])
}

add_network_aracne_edges <- function(graph, aracne_edges, mi_threshold, preserve_duplicates) {
  add <- aracne_edges |>
    dplyr::filter(.data$MI >= mi_threshold, .data$Regulator %in% igraph::V(graph)$name, .data$Target %in% igraph::V(graph)$name)
  if (!nrow(add)) return(graph)
  current <- igraph::as_data_frame(graph, what = "edges") |>
    dplyr::transmute(Regulator = .data$from, Target = .data$to, source = .data$source %||% "regulatory", weight = .data$weight %||% NA_real_)
  add <- add |> dplyr::mutate(weight = .data$MI)
  full <- dplyr::bind_rows(current, add)
  if (!isTRUE(preserve_duplicates)) full <- dplyr::distinct(full)
  igraph::graph_from_data_frame(full |> dplyr::rename(from = "Regulator", to = "Target"), vertices = igraph::as_data_frame(graph, what = "vertices"), directed = FALSE)
}

network_dependencies <- function(networks) {
  unique(c(networks, if ("network_3" %in% networks) "network_1", if ("network_4" %in% networks) "network_2"))
}

write_network_graph <- function(config, graph, network) {
  hotnet <- igraph_to_hotnet(graph)
  save_hotnet_network_files(hotnet, network, out_dir = analysis_paths(config)[["hhotnet_data"]])
  metadata <- igraph::as_data_frame(graph, what = "edges")
  readr::write_tsv(metadata, network_data_file(config, paste0(network, "_edge_metadata.tsv")))
  invalid <- anyDuplicated(hotnet$index$vertex) || any(!hotnet$edge_full$from %in% hotnet$index$vertex) || any(!hotnet$edge_full$to %in% hotnet$index$vertex)
  if (invalid) stop(network, " failed node-ID or endpoint validation", call. = FALSE)
  data.frame(network = network, vertices = igraph::vcount(graph), edges = igraph::ecount(graph),
             components = igraph::components(graph)$no, largest_component = max(igraph::components(graph)$csize),
             invalid_node_or_endpoint = FALSE)
}

read_network_index <- function(file) {
  readr::read_tsv(file, col_names = c("id", "gene"),
                  col_types = readr::cols(id = readr::col_integer(), gene = readr::col_character()),
                  show_col_types = FALSE, progress = FALSE)
}

network_edge_keys <- function(file) {
  edges <- readr::read_tsv(file, show_col_types = FALSE, progress = FALSE)
  required <- c("from", "to", "source")
  if (length(setdiff(required, names(edges)))) {
    stop("Network edge metadata lacks: ", paste(setdiff(required, names(edges)), collapse = ", "), call. = FALSE)
  }
  paste(edges$from, edges$to, edges$source, sep = "\r")
}

#' Compare constructed h-HotNet input networks with a read-only reference
#'
#' Compares semantic node names and edge endpoints/source labels, so the
#' comparison is independent of numeric node IDs. The configured reference is
#' never modified.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for the configured default, one network name, or a
#'   vector of network names.
#' @return A comparison data frame written to `hHotnet/summary`.
#' @export
compare_hhotnet_networks <- function(config, networks = NULL) {
  requested <- resolve_networks(config, networks)
  reference_dir <- config$validation$reference_hhotnet_dir %||% NULL
  if (is.null(reference_dir)) stop("validation.reference_hhotnet_dir is required for network comparison", call. = FALSE)
  reference_dir <- resolve_config_path(reference_dir, config$project_root)
  if (!dir.exists(reference_dir)) stop("Missing reference h-HotNet directory: ", reference_dir, call. = FALSE)
  rows <- lapply(requested, function(network) {
    observed_index_file <- network_data_file(config, paste0(network, "_index_gene.tsv"))
    observed_edge_file <- network_data_file(config, paste0(network, "_edge_list_full.tsv"))
    reference_index_file <- file.path(reference_dir, "data", paste0(network, "_index_gene.tsv"))
    reference_edge_file <- file.path(reference_dir, "data", paste0(network, "_edge_list_full.tsv"))
    files <- c(observed_index_file, observed_edge_file, reference_index_file, reference_edge_file)
    if (any(!file.exists(files))) stop("Missing network comparison input(s): ", paste(files[!file.exists(files)], collapse = ", "), call. = FALSE)
    observed_nodes <- read_network_index(observed_index_file)$gene
    reference_nodes <- read_network_index(reference_index_file)$gene
    observed_edges <- network_edge_keys(observed_edge_file)
    reference_edges <- network_edge_keys(reference_edge_file)
    data.frame(
      network = network,
      observed_nodes = length(observed_nodes), reference_nodes = length(reference_nodes), shared_nodes = length(intersect(observed_nodes, reference_nodes)),
      observed_only_nodes = length(setdiff(observed_nodes, reference_nodes)), reference_only_nodes = length(setdiff(reference_nodes, observed_nodes)),
      observed_edge_rows = length(observed_edges), reference_edge_rows = length(reference_edges), shared_unique_edges = length(intersect(unique(observed_edges), unique(reference_edges))),
      observed_only_unique_edges = length(setdiff(unique(observed_edges), unique(reference_edges))), reference_only_unique_edges = length(setdiff(unique(reference_edges), unique(observed_edges))),
      stringsAsFactors = FALSE
    )
  })
  report <- dplyr::bind_rows(rows)
  output <- network_summary_file(config, "network_reference_comparison.tsv")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(report, output)
  report
}

#' Construct selected standard h-HotNet input networks
#'
#' Produces only the requested member(s) of Networks 1--4. Dependencies are
#' built in memory as needed (Network 3 needs Network 1; Network 4 needs
#' Network 2), but their files are not written unless requested.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for `networks.selected`, a single network, or a
#'   vector of network names.
#' @param dry_run If `TRUE`, return planned products without loading sources.
#' @return Named paths for a dry run or graph statistics for written networks.
#' @export
build_hhotnet_networks <- function(config, networks = NULL, dry_run = TRUE) {
  requested <- resolve_networks(config, networks)
  mhc <- config$regulatory$mhc_exclusion %||% list(enabled = FALSE)
  if (isTRUE(dry_run)) {
    return(stats::setNames(vapply(requested, function(network) network_data_file(config, paste0(network, "_edge_list.tsv")), character(1)), requested))
  }
  require_network_packages()
  paths <- analysis_paths(config)
  dir.create(paths[["hhotnet_data"]], recursive = TRUE, showWarnings = FALSE)
  dir.create(paths[["hhotnet_summary"]], recursive = TRUE, showWarnings = FALSE)
  data <- load_network_sources(config)
  dedup <- deduplicate_network_jeme(data, isTRUE(config$networks$deduplicate_jeme_edges_before_harmonization))
  split <- split_network_multi_promoters(dedup$data, isTRUE(config$networks$split_multi_promoter_edges))
  harmonized <- harmonize_network_promoters(
    split$data,
    config$regulatory$promoter_harmonization %||% list(),
    config
  )
  classified <- classify_network_enhancer_edges(harmonized$data, config)
  ldak <- filter_analysis_network_ldak(classified$data, config)
  cohort <- config$gwas$dataset
  enhancers <- ldak$enhancers |>
    dplyr::filter(.data$cohort == cohort) |>
    dplyr::mutate(node_type = "enhancer", score = score_from_pvalue(.data[[config$scores$pvalue_column]], config$scores$infinite_score_cap))
  promoters <- ldak$promoters |>
    dplyr::filter(.data$cohort == cohort) |>
    dplyr::mutate(node_type = "promoter_gene", score = score_from_pvalue(.data[[config$scores$pvalue_column]], config$scores$infinite_score_cap))
  mhc_result <- exclude_network_mhc(dplyr::bind_rows(enhancers, promoters), mhc)
  nodes <- mhc_result$data |>
    dplyr::mutate(.input_order = dplyr::row_number()) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene), !is.na(.data$score)) |>
    dplyr::group_by(.data$gene) |>
    dplyr::slice_max(.data$score, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup()
  if (identical(config$networks$node_order, "legacy_input")) nodes <- dplyr::arrange(nodes, .data$.input_order) else nodes <- dplyr::arrange(nodes, dplyr::desc(.data$score), .data$gene)
  nodes <- dplyr::select(nodes, -".input_order")
  if (anyDuplicated(nodes$gene) || any(!is.finite(nodes$score)) || any(nodes$score < 0)) stop("Invalid unique node-score table", call. = FALSE)
  jeme <- dplyr::bind_rows(classified$data$jeme)
  jeme_weight <- if ("score" %in% names(jeme)) suppressWarnings(as.numeric(jeme$score)) else rep(NA_real_, nrow(jeme))
  regulatory_all <- dplyr::bind_rows(
    dplyr::transmute(jeme, Regulator = .data$enhancer, Target = .data$promoter, source = "JEME", weight = jeme_weight),
    dplyr::transmute(classified$data$hic$PO, Regulator = .data$Interacting_fragment, Target = .data$Promoter, source = "HiC_PO", weight = NA_real_),
    dplyr::transmute(classified$data$hic$PP, Regulator = .data$Promoter1, Target = .data$Promoter2, source = "HiC_PP", weight = NA_real_)
  )
  preserve_self_loops <- isTRUE(config$networks$preserve_self_loops %||% TRUE)
  aracne_all <- classified$data$aracne |>
    dplyr::mutate(source = "ARACNe") |>
    dplyr::relocate("Regulator", "Target", "source") |>
    dplyr::filter(!is.na(.data$Regulator), !is.na(.data$Target), nzchar(.data$Regulator), nzchar(.data$Target)) |>
    dplyr::distinct() |>
    apply_network_self_loop_policy(preserve_self_loops)
  filtered <- filter_network_regulatory_edges(
    regulatory_all, nodes, config$networks$regulatory_edge_filter,
    preserve_self_loops
  )
  aracne <- aracne_all |> dplyr::filter(.data$Regulator %in% nodes$gene, .data$Target %in% nodes$gene)
  aracne_harmonization_summary <- harmonized$data$aracne_harmonization$edge_summary |>
    dplyr::mutate(
      edges_after_self_loop_policy = nrow(aracne_all),
      edges_after_network_node_filter = nrow(aracne),
      edges_dropped_by_network_node_filter = nrow(aracne_all) - nrow(aracne)
    )
  excluded_ids <- unique(as.character(mhc_result$excluded$gene))
  regulatory_edges_excluded <- sum(
    regulatory_all$Regulator %in% excluded_ids |
      regulatory_all$Target %in% excluded_ids
  )
  aracne_edges_excluded <- sum(
    aracne_all$Regulator %in% excluded_ids |
      aracne_all$Target %in% excluded_ids
  )
  readr::write_tsv(nodes, network_data_file(config, "nodes_all_info.tsv"))
  readr::write_tsv(dplyr::select(nodes, "gene", "score"), network_data_file(config, paste0(config$hhotnet$score_name, ".tsv")), col_names = FALSE)
  readr::write_tsv(filtered$kept, network_data_file(config, "regulatory_edges.tsv"))
  readr::write_tsv(aracne, network_data_file(config, "aracne_edges.tsv"))
  readr::write_tsv(ldak$unmatched_enhancers, network_summary_file(config, "unmatched_enhancers.tsv"))
  readr::write_tsv(ldak$unmatched_promoters, network_summary_file(config, "unmatched_promoters.tsv"))
  readr::write_tsv(dedup$summary, network_summary_file(config, "jeme_deduplication_summary.tsv"))
  readr::write_tsv(split$summary, network_summary_file(config, "multi_promoter_splitting_summary.tsv"))
  readr::write_tsv(harmonized$summary, network_summary_file(config, "promoter_harmonization_summary.tsv"))
  readr::write_tsv(
    harmonized$jeme_mapping,
    network_summary_file(config, "jeme_all_tissue_ensg_harmonization.tsv.gz")
  )
  readr::write_tsv(classified$audit, network_summary_file(config, "regulatory_edge_classification_runtime.tsv.gz"))
  readr::write_tsv(classified$summary, network_summary_file(config, "enhancer_edge_classification_summary.tsv"))
  readr::write_tsv(
    harmonized$data$aracne_harmonization$mapping,
    network_summary_file(config, "aracne_harmonization_mapping.tsv")
  )
  readr::write_tsv(
    harmonized$data$aracne_harmonization$endpoint_summary,
    network_summary_file(config, "aracne_harmonization_endpoints.tsv")
  )
  readr::write_tsv(
    aracne_harmonization_summary,
    network_summary_file(config, "aracne_harmonization_summary.tsv")
  )
  readr::write_tsv(mhc_result$excluded, network_summary_file(config, "mhc_excluded_nodes.tsv"))
  readr::write_tsv(
    data.frame(
      enabled = isTRUE(mhc$enabled),
      genome_build = mhc_result$region$genome_build,
      chr = mhc_result$region$chr,
      start = mhc_result$region$start,
      end = mhc_result$region$end,
      overlap_rule = "inclusive_any_overlap",
      nodes_excluded = nrow(mhc_result$excluded),
      enhancers_excluded = sum(mhc_result$excluded$node_type == "enhancer"),
      promoters_excluded = sum(mhc_result$excluded$node_type == "promoter_gene"),
      other_nodes_excluded = sum(!mhc_result$excluded$node_type %in% c("enhancer", "promoter_gene")),
      regulatory_edges_excluded = regulatory_edges_excluded,
      aracne_edges_excluded = aracne_edges_excluded
    ),
    network_summary_file(config, "mhc_exclusion_summary.tsv")
  )
  base <- make_network_graph(filtered$kept, nodes)
  needed <- network_dependencies(requested)
  graphs <- list()
  if (any(needed %in% c("network_1", "network_3"))) graphs$network_1 <- subset_network_components(base, config$networks$thresholds$small_regulatory_score)
  if (any(needed %in% c("network_2", "network_4"))) graphs$network_2 <- subset_network_components(base, config$networks$thresholds$broad_regulatory_score)
  if ("network_3" %in% needed) graphs$network_3 <- add_network_aracne_edges(graphs$network_1, aracne, config$networks$thresholds$small_plus_aracne_mi, isTRUE(config$networks$preserve_duplicate_edges_during_aracne))
  if ("network_4" %in% needed) graphs$network_4 <- add_network_aracne_edges(graphs$network_2, aracne, config$networks$thresholds$broad_plus_aracne_mi, isTRUE(config$networks$preserve_duplicate_edges_during_aracne))
  stats <- dplyr::bind_rows(lapply(requested, function(network) write_network_graph(config, graphs[[network]], network)))
  readr::write_tsv(stats, network_summary_file(config, "network_file_summary.tsv"))
  stats
}
