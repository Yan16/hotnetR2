# ARACNe endpoint harmonization and provenance reporting.

default_aracne_harmonization <- function(config = NULL) {
  config <- config %||% list()
  list(
    enabled = config$enabled %||% TRUE,
    gtf_file = config$gtf_file %||% ".cache/gencode.v39.GRCh38.genes_nodup.rds",
    alias_file = config$alias_file %||% ".cache/alias_link_nodup.rds",
    map_ensembl = config$map_ensembl %||% TRUE,
    map_aliases = config$map_aliases %||% TRUE,
    ambiguous = config$ambiguous %||% "drop",
    unmapped_ensembl = config$unmapped_ensembl %||% "drop",
    unmatched_symbols = config$unmatched_symbols %||% "keep"
  )
}

normalise_aracne_identifier <- function(x) {
  x <- trimws(as.character(x))
  is_ensg <- grepl("^ENSG[0-9]+\\.[0-9]+$", x, ignore.case = TRUE)
  x[is_ensg] <- sub("\\.[0-9]+$", "", x[is_ensg])
  toupper(x)
}

aracne_mapping_candidates <- function(settings, project_root) {
  rows <- list()
  if (isTRUE(settings$map_ensembl)) {
    gtf_file <- resolve_config_path(settings$gtf_file, project_root)
    if (!file.exists(gtf_file)) {
      stop("Missing ARACNe GENCODE mapping: ", gtf_file, call. = FALSE)
    }
    gtf <- readRDS(gtf_file)
    required <- c("gene_id", "gene_name")
    if (length(setdiff(required, names(gtf)))) {
      stop("ARACNe GENCODE mapping must contain gene_id and gene_name", call. = FALSE)
    }
    rows$gencode_ensg <- data.frame(
      identifier = gtf$gene_id,
      canonical = gtf$gene_name,
      source = "gencode_ensg"
    )
    rows$gencode_symbol <- data.frame(
      identifier = gtf$gene_name,
      canonical = gtf$gene_name,
      source = "gencode_symbol"
    )
  }
  if (isTRUE(settings$map_aliases)) {
    alias_file <- resolve_config_path(settings$alias_file, project_root)
    if (!file.exists(alias_file)) {
      stop("Missing ARACNe alias mapping: ", alias_file, call. = FALSE)
    }
    aliases <- readRDS(alias_file)
    required <- c("ensg", "gene", "Hsym")
    if (length(setdiff(required, names(aliases)))) {
      stop("ARACNe alias mapping must contain ensg, gene, and Hsym", call. = FALSE)
    }
    rows$alias_ensg <- data.frame(
      identifier = aliases$ensg,
      canonical = aliases$gene,
      source = "alias_ensg"
    )
    rows$hgnc_alias <- data.frame(
      identifier = aliases$Hsym,
      canonical = aliases$gene,
      source = "hgnc_alias"
    )
    rows$canonical_symbol <- data.frame(
      identifier = aliases$gene,
      canonical = aliases$gene,
      source = "canonical_symbol"
    )
    if (identical(attr(aliases, "mapping_direction"), "raw_to_approved")) {
      rows$alias_ensg$canonical <- dplyr::coalesce(aliases$Hsym, aliases$gene)
      rows$hgnc_alias <- dplyr::mutate(alias_symbol_pairs(aliases), source = "hgnc_alias")
      rows$canonical_symbol$identifier <- rows$alias_ensg$canonical
      rows$canonical_symbol$canonical <- rows$alias_ensg$canonical
    }
  }
  candidates <- dplyr::bind_rows(rows) |>
    dplyr::mutate(
      identifier = trimws(as.character(.data$identifier)),
      canonical = trimws(as.character(.data$canonical)),
      key = normalise_aracne_identifier(.data$identifier)
    ) |>
    dplyr::filter(
      !is.na(.data$identifier), nzchar(.data$identifier),
      !is.na(.data$canonical), nzchar(.data$canonical),
      !is.na(.data$key), nzchar(.data$key)
    ) |>
    dplyr::distinct(.data$key, .data$canonical, .data$source)

  candidates |>
    dplyr::summarise(
      candidate_count = dplyr::n_distinct(.data$canonical),
      candidates = paste(sort(unique(.data$canonical)), collapse = ";"),
      canonical = if (dplyr::n_distinct(.data$canonical) == 1L) {
        dplyr::first(.data$canonical)
      } else {
        NA_character_
      },
      mapping_source = paste(sort(unique(.data$source)), collapse = ";"),
      .by = "key"
    )
}

map_aracne_identifiers <- function(x, mapping) {
  original <- trimws(as.character(x))
  key <- normalise_aracne_identifier(original)
  index <- match(key, mapping$key)
  candidate_count <- mapping$candidate_count[index]
  mapped <- mapping$canonical[index]
  candidates <- mapping$candidates[index]
  source <- mapping$mapping_source[index]
  is_ensg <- grepl("^ENSG[0-9]+$", key)
  ambiguous <- !is.na(candidate_count) & candidate_count > 1L
  found <- !is.na(mapped) & !ambiguous

  canonical <- mapped
  canonical[!found & !is_ensg] <- original[!found & !is_ensg]
  canonical[ambiguous | (!found & is_ensg)] <- NA_character_
  status <- rep("unmatched_symbol_retained", length(original))
  status[found & normalise_aracne_identifier(mapped) == key] <- "identity"
  status[found & is_ensg] <- "ensembl_mapped"
  status[found & is_ensg & normalise_aracne_identifier(mapped) == key] <-
    "ensembl_identity"
  status[found & !is_ensg & normalise_aracne_identifier(mapped) != key] <- "alias_mapped"
  status[ambiguous] <- "ambiguous_dropped"
  status[!found & !ambiguous & is_ensg] <- "unmapped_ensembl_dropped"
  status[is.na(original) | !nzchar(original)] <- "empty_dropped"
  canonical[is.na(original) | !nzchar(original)] <- NA_character_
  source[!found & !is_ensg] <- "retained_as_supplied"
  source[!found & is_ensg] <- "no_mapping"
  source[ambiguous] <- "conflicting_candidates"

  data.frame(
    original = original,
    canonical = canonical,
    mapping_status = status,
    mapping_source = source,
    mapping_candidates = candidates,
    stringsAsFactors = FALSE
  )
}

#' Harmonize ARACNe regulator and target identifiers
#'
#' Ensembl identifiers are mapped to gene symbols and historical HGNC aliases
#' are mapped to the canonical symbols used by the regulatory networks.
#' Ambiguous and unmapped Ensembl identifiers are dropped conservatively;
#' otherwise-unmatched symbols are retained as supplied. The behavior is
#' enabled by default when `aracne.harmonization` is absent from the YAML.
#'
#' @param aracne ARACNe edge data with `Regulator` and `Target` columns.
#' @param config The `aracne.harmonization` configuration block, or `NULL`.
#' @param project_root Project root used to resolve mapping paths.
#' @return A list containing harmonized `data`, identifier-level `mapping`,
#'   and edge/endpoint summaries.
#' @export
harmonize_aracne_endpoints <- function(aracne, config = NULL, project_root = ".") {
  settings <- default_aracne_harmonization(config)
  required <- c("Regulator", "Target")
  if (length(setdiff(required, names(aracne)))) {
    stop("ARACNe data must contain Regulator and Target", call. = FALSE)
  }
  if (!isTRUE(settings$enabled)) {
    return(list(
      data = aracne,
      mapping = data.frame(),
      endpoint_summary = data.frame(endpoint = "both", mapping_status = "disabled",
                                    identifiers = NA_integer_, occurrences = 2L * nrow(aracne)),
      edge_summary = data.frame(
        harmonization_enabled = FALSE, input_edges = nrow(aracne),
        harmonized_edges = nrow(aracne), dropped_edges = 0L,
        changed_edges = 0L
      )
    ))
  }
  mapping <- aracne_mapping_candidates(settings, project_root)
  regulator <- map_aracne_identifiers(aracne$Regulator, mapping)
  target <- map_aracne_identifiers(aracne$Target, mapping)
  output <- aracne
  output$Regulator_original <- regulator$original
  output$Target_original <- target$original
  output$Regulator <- regulator$canonical
  output$Target <- target$canonical
  output$Regulator_mapping_status <- regulator$mapping_status
  output$Target_mapping_status <- target$mapping_status
  output$Regulator_mapping_source <- regulator$mapping_source
  output$Target_mapping_source <- target$mapping_source
  output$Regulator_mapping_candidates <- regulator$mapping_candidates
  output$Target_mapping_candidates <- target$mapping_candidates
  keep <- !is.na(output$Regulator) & nzchar(output$Regulator) &
    !is.na(output$Target) & nzchar(output$Target)
  changed <- keep & (
    output$Regulator != output$Regulator_original |
      output$Target != output$Target_original
  )

  identifier_map <- dplyr::bind_rows(
    dplyr::mutate(regulator, endpoint = "Regulator"),
    dplyr::mutate(target, endpoint = "Target")
  ) |>
    dplyr::distinct(.data$endpoint, .data$original, .data$canonical,
                    .data$mapping_status, .data$mapping_source,
                    .data$mapping_candidates)
  endpoint_summary <- dplyr::bind_rows(
    dplyr::mutate(regulator, endpoint = "Regulator"),
    dplyr::mutate(target, endpoint = "Target")
  ) |>
    dplyr::summarise(
      identifiers = dplyr::n_distinct(.data$original),
      occurrences = dplyr::n(),
      .by = c("endpoint", "mapping_status")
    ) |>
    dplyr::arrange(.data$endpoint, .data$mapping_status)
  edge_summary <- data.frame(
    harmonization_enabled = TRUE,
    input_edges = nrow(aracne),
    harmonized_edges = sum(keep),
    dropped_edges = sum(!keep),
    changed_edges = sum(changed),
    regulator_identifiers_changed = sum(
      regulator$canonical != regulator$original, na.rm = TRUE
    ),
    target_identifiers_changed = sum(
      target$canonical != target$original, na.rm = TRUE
    )
  )
  list(
    data = output[keep, , drop = FALSE],
    mapping = identifier_map,
    endpoint_summary = endpoint_summary,
    edge_summary = edge_summary
  )
}
