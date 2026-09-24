# LDAK annotation preparation using all tissues from each enabled source.

collapse_annotation_values <- function(x) {
  paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")
}

normalise_ldak_chr <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  x[toupper(x) == "X"] <- "23"
  x[toupper(x) == "Y"] <- "24"
  out <- suppressWarnings(as.integer(x))
  if (anyNA(out)) stop("Annotation contains a non-standard chromosome", call. = FALSE)
  out
}

normalise_genome_build <- function(x) {
  build <- toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
  aliases <- c(HG19 = "GRCH37", GRCH37 = "GRCH37", HG38 = "GRCH38", GRCH38 = "GRCH38")
  out <- unname(aliases[build])
  if (length(out) != 1L || is.na(out)) {
    stop("Genome build must be one of: GRCh37, hg19, GRCh38, hg38", call. = FALSE)
  }
  out
}

normalise_ldak_coordinate <- function(x, field) {
  out <- suppressWarnings(as.numeric(x))
  if (anyNA(out) || any(!is.finite(out)) || any(out != floor(out))) {
    stop("Annotation contains a non-integer or missing ", field, " coordinate", call. = FALSE)
  }
  as.integer(out)
}

ldak_reference_bounds <- function(config) {
  bim <- paste0(resolve_config_path(config$reference$bfile_prefix, config$project_root), ".bim")
  if (!file.exists(bim)) stop("Missing LDAK reference BIM file: ", bim, call. = FALSE)
  bim_data <- readr::read_tsv(
    bim, col_names = FALSE, col_select = c("X1", "X4"),
    col_types = readr::cols(X1 = readr::col_integer(), X4 = readr::col_double()),
    progress = FALSE, show_col_types = FALSE
  )
  dplyr::summarise(
    dplyr::group_by(bim_data, .data$X1),
    min_bp = min(.data$X4, na.rm = TRUE), max_bp = max(.data$X4, na.rm = TRUE),
    reference_variants = dplyr::n(), .groups = "drop"
  ) |>
    dplyr::rename(CHR = .data$X1)
}

validate_ldak_intervals <- function(intervals, annotation_type, reference_bounds) {
  required <- c("name", "CHR", "START", "END")
  missing <- setdiff(required, names(intervals))
  if (length(missing)) {
    stop(annotation_type, " intervals lack column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data <- dplyr::mutate(
    intervals,
    name = trimws(as.character(.data$name)),
    CHR = normalise_ldak_chr(.data$CHR),
    START = normalise_ldak_coordinate(.data$START, "START"),
    END = normalise_ldak_coordinate(.data$END, "END")
  ) |>
    dplyr::relocate("name", "CHR", "START", "END")
  if (any(!nzchar(data$name))) stop(annotation_type, " intervals contain an empty name", call. = FALSE)
  if (any(data$CHR < 1L | data$CHR > 24L)) {
    stop(annotation_type, " intervals must use autosomes or chromosomes X/Y", call. = FALSE)
  }
  if (any(data$START < 0L | data$END <= data$START)) {
    stop(annotation_type, " intervals have invalid START/END coordinates", call. = FALSE)
  }
  if (anyDuplicated(data$name)) {
    stop(annotation_type, " intervals contain duplicate names", call. = FALSE)
  }
  if (anyDuplicated(data[c("name", "CHR", "START", "END")])) {
    stop(annotation_type, " intervals contain duplicate locations", call. = FALSE)
  }
  bounded <- dplyr::left_join(data, reference_bounds, by = "CHR")
  # BIM positions are 1-based points; BED [START,END) covers biological bases
  # START+1 through END. Therefore START == max_bp lies immediately to the
  # right of the final reference position.
  outside <- is.na(bounded$min_bp) | bounded$END < bounded$min_bp | bounded$START >= bounded$max_bp
  report <- data.frame(
    annotation_type = annotation_type,
    metric = c("interval_rows", "unique_names", "chromosomes", "reference_panel_chromosomes",
               "intervals_overlapping_reference_panel", "intervals_fully_outside_reference_panel"),
    value = c(nrow(data), dplyr::n_distinct(data$name), dplyr::n_distinct(data$CHR),
              nrow(reference_bounds), sum(!outside), sum(outside))
  )
  list(data = dplyr::select(data, dplyr::everything()), report = report)
}

annotation_provenance_summary <- function(inputs, annotation_type) {
  dplyr::mutate(inputs, annotation_type = annotation_type) |>
    dplyr::group_by(.data$annotation_type, .data$source, .data$interaction_type) |>
    dplyr::summarise(
      input_rows = dplyr::n(), input_names = dplyr::n_distinct(.data$name),
      input_tissues = dplyr::n_distinct(.data$tissues), .groups = "drop"
    ) |>
    dplyr::arrange(.data$annotation_type, .data$source, .data$interaction_type)
}

require_annotation_packages <- function() {
  packages <- c("dplyr", "readr", "stringr", "tidyr")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("LDAK annotation preparation requires: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

ldak_annotation_config <- function(config) {
  settings <- config$ldak$annotation_tissues %||% list()
  jeme_tissues <- settings$jeme %||% "ALL"
  hic_tissues <- settings$hic %||% "ALL"
  if (!identical(toupper(as.character(jeme_tissues)), "ALL") ||
      (analysis_hic_enabled(config) && !identical(toupper(as.character(hic_tissues)), "ALL"))) {
    stop(
      "The initial prepare_ldak_annotations() implementation supports only ",
      "ldak.annotation_tissues.jeme: ALL and hic: ALL when HiC is enabled", call. = FALSE
    )
  }
  list(jeme_method = settings$jeme_method %||% "lasso")
}

ncbi_feature_table_path <- function(config) {
  path <- config$reference$ncbi_feature_table
  if (is.null(path)) {
    stop(
      "reference.ncbi_feature_table must name an existing feature table before ",
      "promoter annotations can be created. Automatic downloads are disabled ",
      "to preserve the analysis1 development boundary.",
      call. = FALSE
    )
  }
  path <- resolve_config_path(path, config$project_root)
  if (!file.exists(path)) stop("Missing reference.ncbi_feature_table: ", path, call. = FALSE)
  path
}

promoter_alias_map_path <- function(config) {
  path <- config$reference$promoter_alias_file %||% NULL
  if (is.null(path)) return(NULL)
  path <- resolve_config_path(path, config$project_root)
  if (!file.exists(path)) stop("Missing reference.promoter_alias_file: ", path, call. = FALSE)
  path
}

harmonize_promoter_aliases <- function(promoter_sources, config) {
  alias_path <- promoter_alias_map_path(config)
  if (is.null(alias_path)) {
    return(dplyr::mutate(
      promoter_sources, input_symbols = .data$name,
      alias_harmonization = "identity_no_alias_file"
    ))
  }
  columns <- config$reference$promoter_alias_columns %||% list(alias = "alias", gene_symbol = "gene_symbol")
  for (field in c("alias", "gene_symbol")) {
    assert_scalar(columns[[field]], paste0("reference.promoter_alias_columns.", field), "character")
  }
  aliases <- readr::read_tsv(alias_path, show_col_types = FALSE, progress = FALSE)
  required <- unname(unlist(columns[c("alias", "gene_symbol")]))
  missing <- setdiff(required, names(aliases))
  if (length(missing)) {
    stop("Promoter alias file lacks configured column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  aliases <- dplyr::transmute(
    aliases,
    alias = trimws(as.character(.data[[columns$alias]])),
    gene_symbol = trimws(as.character(.data[[columns$gene_symbol]]))
  ) |>
    dplyr::filter(nzchar(.data$alias), nzchar(.data$gene_symbol)) |>
    dplyr::distinct()
  conflicting <- dplyr::count(aliases, .data$alias, name = "n") |>
    dplyr::filter(.data$n > 1L)
  if (nrow(conflicting)) {
    stop("Promoter alias file maps at least one alias to multiple gene symbols", call. = FALSE)
  }
  aliases <- dplyr::bind_rows(
    aliases,
    dplyr::transmute(aliases, alias = .data$gene_symbol, gene_symbol = .data$gene_symbol)
  ) |>
    dplyr::distinct(.data$alias, .keep_all = TRUE)

  dplyr::left_join(promoter_sources, aliases, by = c("name" = "alias")) |>
    dplyr::mutate(
      input_symbols = .data$name,
      alias_harmonization = ifelse(is.na(.data$gene_symbol), "unmapped_alias", "mapped_alias"),
      name = dplyr::coalesce(.data$gene_symbol, .data$name)
    ) |>
    dplyr::select(-.data$gene_symbol) |>
    dplyr::group_by(.data$name) |>
    dplyr::summarise(
      source = collapse_annotation_values(.data$source),
      interaction_type = collapse_annotation_values(.data$interaction_type),
      tissues = collapse_annotation_values(unlist(stringr::str_split(.data$tissues, "[,;]\\s*"))),
      input_symbols = collapse_annotation_values(.data$input_symbols),
      alias_harmonization = collapse_annotation_values(.data$alias_harmonization),
      .groups = "drop"
    )
}

convert_promoter_details_to_tss <- function(promoter_details) {
  required <- c("START", "END", "strand")
  missing <- setdiff(required, names(promoter_details))
  if (length(missing)) {
    stop("Promoter details lack: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invalid_strand <- !promoter_details$strand %in% c("+", "-")
  if (any(invalid_strand)) {
    stop(
      "Cannot construct strand-aware TSS for ", sum(invalid_strand),
      " promoter row(s) with unsupported strand", call. = FALSE
    )
  }
  promoter_details |>
    dplyr::mutate(
      gene_START = .data$START,
      gene_END = .data$END,
      # NCBI source bounds are 1-based inclusive. Store the biological TSS as
      # a 1-based base position and write its canonical one-base BED interval.
      TSS = dplyr::if_else(.data$strand == "+", .data$START, .data$END),
      START = .data$TSS - 1L,
      END = .data$TSS,
      promoter_interval_mode = "strand_aware_tss"
    ) |>
    dplyr::relocate(
      "name", "CHR", "START", "END", "strand", "TSS",
      "gene_START", "gene_END", "promoter_interval_mode"
    ) |>
    dplyr::arrange(.data$CHR, .data$START, .data$END, .data$name)
}

convert_promoter_details_to_gene_body <- function(promoter_details) {
  required <- c("START", "END")
  missing <- setdiff(required, names(promoter_details))
  if (length(missing)) {
    stop("Promoter details lack: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  promoter_details |>
    dplyr::mutate(
      gene_START = .data$START,
      gene_END = .data$END,
      # Convert NCBI's 1-based inclusive gene bounds to BED exactly once.
      START = .data$gene_START - 1L,
      END = .data$gene_END,
      promoter_interval_mode = "gene_body"
    ) |>
    dplyr::relocate(
      "name", "CHR", "START", "END", "strand",
      "gene_START", "gene_END", "promoter_interval_mode"
    ) |>
    dplyr::arrange(.data$CHR, .data$START, .data$END, .data$name)
}

annotation_output_file <- function(config, name) {
  path <- file.path(analysis_paths(config)[["ldak_annotations"]], name)
  root <- analysis_paths(config)[["root"]]
  if (!startsWith(normalizePath(path, mustWork = FALSE), paste0(root, .Platform$file.sep))) {
    stop("Refusing annotation output outside analysis directory", call. = FALSE)
  }
  path
}

#' Validate existing LDAK annotation products
#'
#' Rechecks `.loc`/detail-table mappings, interval order, uniqueness, and
#' reference-panel overlap without retrieving JEME or Hi-C source data. This is
#' useful for validating a completed annotation generation on a workstation
#' with less memory than the full all-tissue extraction requires.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @return Paths to the refreshed validation and provenance reports.
#' @export
validate_ldak_annotations <- function(config) {
  require_annotation_packages()
  specs <- list(
    enhancer = c(loc = annotation_output_file(config, "enhancers.loc"),
                 details = annotation_output_file(config, "enhancers.details.tsv.gz")),
    promoter = c(loc = annotation_output_file(config, "promoters.loc"),
                 details = annotation_output_file(config, "promoters.details.tsv.gz"))
  )
  if (any(!file.exists(unlist(specs)))) {
    missing <- unlist(specs)[!file.exists(unlist(specs))]
    stop("Missing LDAK annotation product(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  bounds <- ldak_reference_bounds(config)
  validation <- list()
  provenance <- list()
  for (type in names(specs)) {
    detail <- readr::read_tsv(specs[[type]][["details"]], show_col_types = FALSE, progress = FALSE)
    checked <- validate_ldak_intervals(detail, type, bounds)
    loc <- readr::read_tsv(
      specs[[type]][["loc"]], col_names = c("name", "CHR", "START", "END"),
      col_types = readr::cols(name = readr::col_character(), CHR = readr::col_integer(),
                              START = readr::col_integer(), END = readr::col_integer()),
      show_col_types = FALSE, progress = FALSE
    )
    loc_matches_detail <- identical(
      as.data.frame(loc),
      as.data.frame(dplyr::select(checked$data, "name", "CHR", "START", "END"))
    )
    expected_order <- order(checked$data$CHR, checked$data$START, checked$data$END, checked$data$name)
    sorted <- identical(expected_order, seq_len(nrow(checked$data)))
    if (!loc_matches_detail || !sorted) {
      stop(type, " .loc file does not exactly match a sorted detail table", call. = FALSE)
    }
    validation[[type]] <- dplyr::bind_rows(
      checked$report,
      data.frame(annotation_type = type, metric = c("loc_matches_detail", "loc_sorted"), value = c(1L, 1L))
    )
    provenance[[type]] <- checked$data |>
      dplyr::group_by(.data$source, .data$interaction_type) |>
      dplyr::summarise(final_loc_rows = dplyr::n(), final_loc_names = dplyr::n_distinct(.data$name), .groups = "drop") |>
      dplyr::mutate(annotation_type = type, provenance_level = "final_deduplicated") |>
      dplyr::relocate(.data$annotation_type, .data$provenance_level)
  }
  outputs <- c(
    interval_validation = annotation_output_file(config, "interval_validation.tsv"),
    provenance_summary = annotation_output_file(config, "annotation_provenance_summary.tsv")
  )
  readr::write_tsv(dplyr::bind_rows(validation), outputs[["interval_validation"]])
  readr::write_tsv(dplyr::bind_rows(provenance), outputs[["provenance_summary"]])
  outputs
}

#' Create all-tissue LDAK enhancer and promoter annotations
#'
#' Uses `hotnetR` to retrieve the same all-JEME/all-Hi-C inputs as the BEEA
#' reference workflow. All outputs are written below `analysis1/LDAK/`.
#' A pre-existing NCBI feature table must be configured; this function never
#' downloads reference data or writes outside the analysis directory.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param dry_run If `TRUE`, validate configuration and return planned outputs
#'   without retrieving data or writing files.
#' @return Named paths and (when run) row counts.
#' @export
prepare_ldak_annotations <- function(config, dry_run = TRUE) {
  settings <- ldak_annotation_config(config)
  outputs <- c(
    enhancer_loc = annotation_output_file(config, "enhancers.loc"),
    enhancer_details = annotation_output_file(config, "enhancers.details.tsv.gz"),
    promoter_loc = annotation_output_file(config, "promoters.loc"),
    promoter_details = annotation_output_file(config, "promoters.details.tsv.gz"),
    jeme_ensg_harmonization = annotation_output_file(config, "jeme_ensg_harmonization.tsv.gz"),
    promoter_harmonization = annotation_output_file(config, "promoter_harmonization_summary.tsv"),
    interval_validation = annotation_output_file(config, "interval_validation.tsv"),
    provenance_summary = annotation_output_file(config, "annotation_provenance_summary.tsv"),
    summary = annotation_output_file(config, "annotation_summary.tsv")
  )
  if (isTRUE(dry_run)) return(outputs)

  require_annotation_packages()
  ncbi_path <- ncbi_feature_table_path(config)
  genome_build <- normalise_genome_build(config$reference$genome_build)
  annotation_genome_build <- normalise_genome_build(config$reference$annotation_genome_build)
  dir.create(analysis_paths(config)[["ldak_annotations"]], recursive = TRUE, showWarnings = FALSE)

  promoter_harmonization <- config$regulatory$promoter_harmonization %||% list()
  if (!isTRUE(promoter_harmonization$enabled %||% FALSE)) {
    stop("All-tissue promoter harmonization must be enabled for v2_analysis2", call. = FALSE)
  }
  gtf_path <- resolve_config_path(
    promoter_harmonization$gtf_file, config$project_root
  )
  hic_alias_path <- if (analysis_hic_enabled(config)) resolve_config_path(
    promoter_harmonization$hic_alias_file, config$project_root
  ) else NULL
  if (!file.exists(gtf_path)) stop("Missing JEME GENCODE mapping: ", gtf_path, call. = FALSE)
  if (!is.null(hic_alias_path) && !file.exists(hic_alias_path)) stop("Missing HiC alias mapping: ", hic_alias_path, call. = FALSE)

  jeme <- get_jeme(method = settings$jeme_method, simplified = TRUE, cache_dir = config$cache_dir)
  jeme_mapping <- build_all_tissue_jeme_ensg_map(jeme, readRDS(gtf_path))
  jeme <- apply_all_tissue_jeme_ensg_map(jeme, jeme_mapping)
  jeme_changed <- sum(jeme$promoter != jeme$promoter_original, na.rm = TRUE)
  jeme_nodes <- extract_jeme_nodes(jeme)
  hic <- load_analysis_hic(config, all_tissues = TRUE)
  hic_harmonized <- harmonize_hic_gene_symbols(hic, hic_alias_path)
  hic <- hic_harmonized$data
  hic_nodes <- extract_hic_nodes(hic)
  input_counts <- list(
    jeme_tissues = dplyr::n_distinct(jeme$desc1),
    hic_po_tissues = dplyr::n_distinct(hic$PO$Tissue_type),
    hic_pp_tissues = dplyr::n_distinct(hic$PP$Tissue_type),
    jeme_interactions = nrow(jeme), hic_po_interactions = nrow(hic$PO),
    hic_pp_interactions = nrow(hic$PP)
  )

  jeme_enh <- dplyr::transmute(
    jeme_nodes$enhancers,
    name = .data$enhancer, CHR = normalise_ldak_chr(.data$CHR),
    START = normalise_ldak_coordinate(.data$START, "START"), END = normalise_ldak_coordinate(.data$END, "END"),
    source = "JEME", interaction_type = "JEME", tissues = .data$tissues
  )
  hic_enh <- dplyr::transmute(
    hic_nodes$enhancers,
    name = .data$Interacting_fragment, CHR = normalise_ldak_chr(.data$CHR),
    START = normalise_ldak_coordinate(.data$START, "START"), END = normalise_ldak_coordinate(.data$END, "END"),
    source = "HiC", interaction_type = "HiC_PO", tissues = .data$tissues
  )
  enhancer_provenance <- dplyr::bind_rows(
    annotation_provenance_summary(jeme_enh, "enhancer"),
    annotation_provenance_summary(hic_enh, "enhancer")
  )
  enhancer_details <- dplyr::bind_rows(jeme_enh, hic_enh) |>
    dplyr::group_by(.data$name, .data$CHR, .data$START, .data$END) |>
    dplyr::summarise(
      source = collapse_annotation_values(.data$source),
      interaction_type = collapse_annotation_values(.data$interaction_type),
      tissues = collapse_annotation_values(unlist(stringr::str_split(.data$tissues, "[,;]\\s*"))),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$CHR, .data$START, .data$END, .data$name)

  jeme_prom <- dplyr::transmute(
    jeme_nodes$promoters, name = .data$promoter, source = "JEME",
    interaction_type = "JEME", tissues = .data$tissues
  )
  hic_po_prom <- tidyr::separate_rows(hic$PO, "Promoter", sep = ";") |>
    dplyr::transmute(
      name = stringr::str_trim(.data$Promoter), source = "HiC",
      interaction_type = "HiC_PO", tissues = .data$Tissue_type
    )
  hic_pp_prom <- dplyr::bind_rows(
    dplyr::transmute(hic$PP, name = .data$Promoter1, tissues = .data$Tissue_type),
    dplyr::transmute(hic$PP, name = .data$Promoter2, tissues = .data$Tissue_type)
  ) |>
    tidyr::separate_rows("name", sep = ";") |>
    dplyr::transmute(
      name = stringr::str_trim(.data$name), source = "HiC",
      interaction_type = "HiC_PP", tissues = .data$tissues
    )
  promoter_inputs <- dplyr::bind_rows(jeme_prom, hic_po_prom, hic_pp_prom) |>
    dplyr::filter(!is.na(.data$name), nzchar(.data$name))
  promoter_provenance <- annotation_provenance_summary(promoter_inputs, "promoter")
  rm(jeme_enh, hic_enh, jeme_nodes, hic_nodes,
     jeme_prom, hic_po_prom, hic_pp_prom, jeme, hic)
  gc(verbose = FALSE)
  promoter_sources <- promoter_inputs |>
    dplyr::filter(!is.na(.data$name), nzchar(.data$name)) |>
    dplyr::group_by(.data$name) |>
    dplyr::summarise(
      source = collapse_annotation_values(.data$source),
      interaction_type = collapse_annotation_values(.data$interaction_type),
      tissues = collapse_annotation_values(unlist(stringr::str_split(.data$tissues, "[,;]\\s*"))),
      .groups = "drop"
    )

  ncbi <- NCBI_feature_table(ncbi_path)
  promoter_sources <- harmonize_promoter_aliases(promoter_sources, config)
  promoter_details <- dplyr::inner_join(
    promoter_sources, ncbi, by = c("name" = "gene_symbol")
  ) |>
    dplyr::select("name", "CHR", "START", "END", "strand", "gene_name", "source",
                  "interaction_type", "tissues", "input_symbols", "alias_harmonization",
                  dplyr::everything()) |>
    dplyr::arrange(.data$CHR, .data$START, .data$END, .data$name) |>
    dplyr::distinct(.data$name, .keep_all = TRUE)

  promoter_interval_mode <- config$ldak$promoter_interval_mode %||% "gene_body"
  if (identical(promoter_interval_mode, "strand_aware_tss")) {
    promoter_details <- convert_promoter_details_to_tss(promoter_details)
  } else {
    promoter_details <- convert_promoter_details_to_gene_body(promoter_details)
  }

  reference_bounds <- ldak_reference_bounds(config)
  enhancer_validation <- validate_ldak_intervals(enhancer_details, "enhancer", reference_bounds)
  promoter_validation <- validate_ldak_intervals(promoter_details, "promoter", reference_bounds)
  enhancer_details <- enhancer_validation$data
  promoter_details <- promoter_validation$data

  readr::write_tsv(dplyr::select(enhancer_details, "name", "CHR", "START", "END"),
                    outputs[["enhancer_loc"]], col_names = FALSE)
  readr::write_tsv(enhancer_details, outputs[["enhancer_details"]])
  readr::write_tsv(enhancer_details, annotation_output_file(config, "enhancers.all.details.tsv.gz"))
  readr::write_tsv(dplyr::select(promoter_details, "name", "CHR", "START", "END"),
                    outputs[["promoter_loc"]], col_names = FALSE)
  readr::write_tsv(promoter_details, outputs[["promoter_details"]])
  readr::write_tsv(jeme_mapping, outputs[["jeme_ensg_harmonization"]])
  readr::write_tsv(
    data.frame(
      source = c("JEME", "HiC_PO", "HiC_PP"),
      method = c(
        "all_tissues_ENSG_map", "gene_symbol_alias_mapping",
        "gene_symbol_alias_mapping"
      ),
      input_rows = c(
        input_counts$jeme_interactions,
        hic_harmonized$summary$hic_po_rows,
        hic_harmonized$summary$hic_pp_rows
      ),
      changed_rows = c(
        jeme_changed,
        hic_harmonized$summary$hic_po_changed,
        hic_harmonized$summary$hic_pp_changed
      )
    ),
    outputs[["promoter_harmonization"]]
  )
  readr::write_tsv(dplyr::bind_rows(enhancer_validation$report, promoter_validation$report),
                    outputs[["interval_validation"]])
  readr::write_tsv(dplyr::bind_rows(enhancer_provenance, promoter_provenance),
                    outputs[["provenance_summary"]])

  summary <- data.frame(
    item = c("jeme_tissues", "hic_po_tissues", "hic_pp_tissues", "jeme_interactions",
             "hic_po_interactions", "hic_pp_interactions", "enhancer_nodes",
             "promoter_symbols_requested", "promoter_nodes_mapped_ncbi",
             "promoter_symbols_unmapped_ncbi", "reference_genome_build",
             "annotation_genome_build", "ldak_enhancer_flank_bp",
             "ldak_promoter_flank_bp", "jeme_ensg_ids_in_all_tissue_map",
             "jeme_rows_changed_by_ensg", "hic_po_rows_changed_by_symbol",
             "hic_pp_rows_changed_by_symbol", "promoter_interval_mode",
             "plus_strand_tss", "minus_strand_tss"),
    value = c(
      input_counts$jeme_tissues, input_counts$hic_po_tissues, input_counts$hic_pp_tissues,
      input_counts$jeme_interactions, input_counts$hic_po_interactions, input_counts$hic_pp_interactions,
      nrow(enhancer_details), nrow(promoter_sources), nrow(promoter_details),
      nrow(dplyr::anti_join(promoter_sources, ncbi, by = c("name" = "gene_symbol"))),
      genome_build, annotation_genome_build,
      as.character(ldak_flank_bp(config, "enhancer")),
      as.character(ldak_flank_bp(config, "promoter")),
      nrow(jeme_mapping), jeme_changed,
      hic_harmonized$summary$hic_po_changed,
      hic_harmonized$summary$hic_pp_changed,
      promoter_interval_mode,
      sum(promoter_details$strand == "+"),
      sum(promoter_details$strand == "-")
    )
  )
  if (identical(promoter_interval_mode, "gene_body")) {
    summary <- dplyr::filter(summary, !.data$item %in% c("promoter_interval_mode", "plus_strand_tss", "minus_strand_tss"))
  }
  readr::write_tsv(summary, outputs[["summary"]])
  c(outputs, enhancer_nodes = nrow(enhancer_details), promoter_nodes = nrow(promoter_details))
}
