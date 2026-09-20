# Configuration and path helpers for the analysis1-local development package.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\]|~)", path)
}

assert_scalar <- function(x, field, type = NULL, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(invisible(TRUE))
    stop("Missing required configuration field: ", field, call. = FALSE)
  }
  if (length(x) != 1L || is.na(x)) {
    stop("Configuration field must be one non-missing value: ", field, call. = FALSE)
  }
  if (!is.null(type) && !inherits(x, type)) {
    stop("Configuration field ", field, " must be a ", type, call. = FALSE)
  }
  invisible(TRUE)
}

assert_flag <- function(x, field, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(TRUE))
  assert_scalar(x, field, "logical")
}

assert_positive_integer <- function(x, field, allow_zero = FALSE) {
  assert_scalar(x, field)
  if (!is.numeric(x) || !is.finite(x) || x != as.integer(x) ||
      (!allow_zero && x < 1) || (allow_zero && x < 0)) {
    qualifier <- if (allow_zero) "a non-negative integer" else "a positive integer"
    stop("Configuration field ", field, " must be ", qualifier, call. = FALSE)
  }
  invisible(TRUE)
}

ldak_flank_bp <- function(config, node_type) {
  if (!node_type %in% c("enhancer", "promoter")) {
    stop("Unsupported LDAK node type: ", node_type, call. = FALSE)
  }
  specific <- config$ldak[[paste0(node_type, "_flank_bp")]]
  value <- specific %||% config$ldak$flank_bp
  if (is.null(value)) {
    stop(
      "Missing LDAK flank for ", node_type, ": configure ldak.",
      node_type, "_flank_bp or the legacy ldak.flank_bp fallback",
      call. = FALSE
    )
  }
  as.integer(value)
}

ldak_result_subdir <- function(config, node_type) {
  if (!node_type %in% c("enhancer", "promoter")) {
    stop("Unsupported LDAK node type: ", node_type, call. = FALSE)
  }
  config$ldak$result_dirs[[node_type]] %||% node_type
}

normalise_network_names <- function(networks, field) {
  if (is.null(networks) || !length(networks)) {
    stop("Configuration field ", field, " must name at least one network", call. = FALSE)
  }
  networks <- trimws(as.character(unlist(networks, use.names = FALSE)))
  if (anyNA(networks) || any(!nzchar(networks))) {
    stop("Configuration field ", field, " contains an empty network name", call. = FALSE)
  }
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", networks))) {
    stop("Network names in ", field, " must be safe identifiers", call. = FALSE)
  }
  if (anyDuplicated(networks)) {
    stop("Configuration field ", field, " contains duplicate network name(s)", call. = FALSE)
  }
  networks
}

resolve_config_path <- function(path, base_dir) {
  assert_scalar(path, "path", "character")
  if (is_absolute_path(path)) {
    return(normalizePath(path, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, path), mustWork = FALSE)
}

mhc_exclusion_region <- function(config = list()) {
  region <- list(
    genome_build = "hg19",
    chr = 6L,
    start = 25000000L,
    end = 34000000L
  )
  if (!is.null(config$genome_build) &&
      !identical(normalise_genome_build(config$genome_build),
                 normalise_genome_build(region$genome_build))) {
    stop(
      "regulatory.mhc_exclusion.genome_build must be hg19/GRCh37",
      call. = FALSE
    )
  }
  for (field in c("chr", "start", "end")) {
    if (!is.null(config[[field]]) &&
        !identical(as.numeric(config[[field]]), as.numeric(region[[field]]))) {
      stop(
        "regulatory.mhc_exclusion.", field, " is fixed at ",
        format(region[[field]], scientific = FALSE, trim = TRUE),
        " for the conservative hg19 MHC exclusion",
        call. = FALSE
      )
    }
  }
  region
}

validate_analysis_config <- function(config) {
  mode <- config$hhotnet$analysis_mode %||% "permutation"
  if (length(mode) != 1L || !mode %in% c("permutation", "manual_delta")) {
    stop("hhotnet.analysis_mode must be permutation or manual_delta", call. = FALSE)
  }
  if (identical(mode, "manual_delta") && !identical(config$hhotnet$execution_mode, "local_python")) {
    stop("manual_delta currently requires local_python execution", call. = FALSE)
  }
  analysis_hic_edge_types(config)
  validate_hic_jeme_overlap(config)
  required_sections <- c(
    "schema_version", "analysis", "project", "gwas", "reference",
    "regulatory", "aracne", "ldak", "scores", "networks", "hhotnet"
  )
  missing <- setdiff(required_sections, names(config))
  if (length(missing)) {
    stop("Configuration missing section(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(config$schema_version) || length(config$schema_version) != 1L ||
      is.na(config$schema_version) || !identical(as.integer(config$schema_version), 1L)) {
    stop("Unsupported schema_version; expected 1", call. = FALSE)
  }
  assert_scalar(config$analysis$name, "analysis.name", "character")
  assert_scalar(config$project$root, "project.root", "character")
  assert_scalar(config$gwas$dataset, "gwas.dataset", "character")
  assert_positive_integer(config$gwas$sample_size, "gwas.sample_size")
  for (field in c("summary_file", "pvalue_file", "extract_file")) {
    assert_scalar(config$gwas[[field]], paste0("gwas.", field), "character")
  }
  for (field in c("ldak_exe", "bfile_prefix")) {
    assert_scalar(config$reference[[field]], paste0("reference.", field), "character")
  }
  for (field in c("genome_build", "annotation_genome_build")) {
    assert_scalar(config$reference[[field]], paste0("reference.", field), "character")
  }
  if (!identical(normalise_genome_build(config$reference$genome_build),
                 normalise_genome_build(config$reference$annotation_genome_build))) {
    stop(
      "reference.genome_build and reference.annotation_genome_build must describe the same build",
      call. = FALSE
    )
  }
  mhc <- config$regulatory$mhc_exclusion %||% list(enabled = FALSE)
  promoter_harmonization <- config$regulatory$promoter_harmonization %||%
    list(enabled = FALSE)
  assert_flag(
    promoter_harmonization$enabled %||% FALSE,
    "regulatory.promoter_harmonization.enabled"
  )
  if (isTRUE(promoter_harmonization$enabled)) {
    for (field in c("jeme", "hic")) {
      assert_flag(
        promoter_harmonization[[field]] %||% TRUE,
        paste0("regulatory.promoter_harmonization.", field)
      )
    }
    mapping_fields <- c("gtf_file", if (analysis_hic_enabled(config)) "hic_alias_file")
    for (field in mapping_fields) {
      assert_scalar(
        promoter_harmonization[[field]],
        paste0("regulatory.promoter_harmonization.", field), "character"
      )
    }
    if (!identical(
      promoter_harmonization$jeme_strategy %||% "all_tissues_ensg_map",
      "all_tissues_ensg_map"
    )) {
      stop("JEME promoter harmonization must use all_tissues_ensg_map", call. = FALSE)
    }
    if (!identical(
      promoter_harmonization$hic_strategy %||% "gene_symbol_alias",
      "gene_symbol_alias"
    )) {
      stop("HiC promoter harmonization must use gene_symbol_alias", call. = FALSE)
    }
  }
  assert_flag(mhc$enabled %||% FALSE, "regulatory.mhc_exclusion.enabled")
  mhc_exclusion_region(mhc)
  if (isTRUE(mhc$enabled)) {
    if (!identical(normalise_genome_build(config$reference$annotation_genome_build),
                   "GRCH37")) {
      stop(
        "MHC exclusion requires hg19/GRCh37 annotation coordinates",
        call. = FALSE
      )
    }
    legacy_fields <- intersect(c("enhancers", "promoters"), names(mhc))
    if (length(legacy_fields)) {
      stop(
        "MHC exclusion applies to all node types; remove legacy field(s): ",
        paste0("regulatory.mhc_exclusion.", legacy_fields, collapse = ", "),
        call. = FALSE
      )
    }
  }
  enhancer_classification <- config$enhancer_classification %||% list(enabled = FALSE)
  assert_flag(
    enhancer_classification$enabled %||% FALSE,
    "enhancer_classification.enabled"
  )
  if (isTRUE(enhancer_classification$enabled)) {
    if (any(!enhancer_drop_sources(config) %in% c("JEME", "HiC_PO"))) {
      stop("drop_non_class1_edge_sources must contain only JEME or HiC_PO", call. = FALSE)
    }
    assert_scalar(
      enhancer_classification$source_analysis,
      "enhancer_classification.source_analysis", "character"
    )
    for (field in c("enhancer_flank_bp", "promoter_flank_bp")) {
      assert_positive_integer(
        enhancer_classification[[field]],
        paste0("enhancer_classification.", field), allow_zero = TRUE
      )
    }
    for (field in c("class1_label", "class2_label", "class3_label")) {
      assert_scalar(
        enhancer_classification[[field]],
        paste0("enhancer_classification.", field), "character"
      )
      if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", enhancer_classification[[field]])) {
        stop("Enhancer class labels must be safe identifiers", call. = FALSE)
      }
    }
    class_labels <- unlist(enhancer_classification[c(
      "class1_label", "class2_label", "class3_label"
    )], use.names = FALSE)
    if (anyDuplicated(class_labels)) {
      stop("Enhancer class labels must differ", call. = FALSE)
    }
    if (!enhancer_classification$promoter_window %in% c("strand_aware_tss", "gene_body")) {
      stop("enhancer_classification.promoter_window must be gene_body or strand_aware_tss",
           call. = FALSE)
    }
  }
  aracne_enabled <- config$aracne$enabled %||% TRUE
  assert_flag(aracne_enabled, "aracne.enabled")
  if (isTRUE(aracne_enabled)) {
    assert_scalar(config$aracne$file, "aracne.file", "character")
    aracne_harmonization <- default_aracne_harmonization(
      config$aracne$harmonization %||% NULL
    )
    for (field in c("enabled", "map_ensembl", "map_aliases")) {
      assert_flag(
        aracne_harmonization[[field]],
        paste0("aracne.harmonization.", field)
      )
    }
    for (field in c("gtf_file", "alias_file")) {
      assert_scalar(
        aracne_harmonization[[field]],
        paste0("aracne.harmonization.", field),
        "character"
      )
    }
    if (!identical(aracne_harmonization$ambiguous, "drop")) {
      stop("aracne.harmonization.ambiguous currently supports only drop",
           call. = FALSE)
    }
    if (!identical(aracne_harmonization$unmapped_ensembl, "drop")) {
      stop("aracne.harmonization.unmapped_ensembl currently supports only drop",
           call. = FALSE)
    }
    if (!identical(aracne_harmonization$unmatched_symbols, "keep")) {
      stop("aracne.harmonization.unmatched_symbols currently supports only keep",
           call. = FALSE)
    }
  }
  if (!is.null(config$ldak$flank_bp)) {
    assert_positive_integer(config$ldak$flank_bp, "ldak.flank_bp", allow_zero = TRUE)
  }
  for (node_type in c("enhancer", "promoter")) {
    specific_field <- paste0(node_type, "_flank_bp")
    if (!is.null(config$ldak[[specific_field]])) {
      assert_positive_integer(
        config$ldak[[specific_field]], paste0("ldak.", specific_field),
        allow_zero = TRUE
      )
    }
    assert_positive_integer(
      ldak_flank_bp(config, node_type), paste0("effective ldak.", specific_field),
      allow_zero = TRUE
    )
    result_subdir <- ldak_result_subdir(config, node_type)
    assert_scalar(result_subdir, paste0("ldak.result_dirs.", node_type), "character")
    if (!grepl("^[A-Za-z][A-Za-z0-9_.-]*$", result_subdir)) {
      stop(
        "ldak.result_dirs.", node_type,
        " must be a safe single directory name", call. = FALSE
      )
    }
  }
  promoter_interval_mode <- config$ldak$promoter_interval_mode %||% "gene_body"
  if (!promoter_interval_mode %in% c("gene_body", "strand_aware_tss")) {
    stop(
      "ldak.promoter_interval_mode must be gene_body or strand_aware_tss",
      call. = FALSE
    )
  }
  assert_positive_integer(config$ldak$permutations, "ldak.permutations", allow_zero = TRUE)
  assert_positive_integer(config$ldak$threads, "ldak.threads")
  assert_flag(config$ldak$ignore_weights, "ldak.ignore_weights")
  assert_flag(config$ldak$allow_ambiguous, "ldak.allow_ambiguous")

  available <- normalise_network_names(config$networks$available %||% config$networks$selected,
                                        "networks.available")
  selected <- normalise_network_names(config$networks$selected, "networks.selected")
  if (any(!selected %in% available)) {
    stop("Every networks.selected value must occur in networks.available", call. = FALSE)
  }
  allowed_filter <- c("current", "legacy_filter_network_ldak")
  if (!config$networks$regulatory_edge_filter %in% allowed_filter) {
    stop("networks.regulatory_edge_filter must be one of: ",
         paste(allowed_filter, collapse = ", "), call. = FALSE)
  }
  allowed_order <- c("score_descending", "legacy_input")
  if (!config$networks$node_order %in% allowed_order) {
    stop("networks.node_order must be one of: ", paste(allowed_order, collapse = ", "), call. = FALSE)
  }
  assert_flag(config$networks$preserve_self_loops %||% TRUE,
              "networks.preserve_self_loops")
  for (field in c("deduplicate_jeme_edges_before_harmonization", "split_multi_promoter_edges")) {
    assert_flag(config$networks[[field]], paste0("networks.", field))
  }
  if (isTRUE(aracne_enabled)) {
    assert_flag(config$networks$preserve_duplicate_edges_during_aracne,
                "networks.preserve_duplicate_edges_during_aracne")
  }
  threshold_fields <- c("small_regulatory_score", "broad_regulatory_score")
  if (isTRUE(aracne_enabled)) {
    threshold_fields <- c(threshold_fields, "small_plus_aracne_mi", "broad_plus_aracne_mi")
  }
  for (field in threshold_fields) {
    value <- config$networks$thresholds[[field]]
    assert_scalar(value, paste0("networks.thresholds.", field))
    if (!is.numeric(value) || !is.finite(value) || value < 0) {
      stop("Network threshold must be a non-negative finite number: ", field, call. = FALSE)
    }
  }
  assert_scalar(config$scores$pvalue_column, "scores.pvalue_column", "character")
  if (!config$scores$transform %in% c("neg_log10")) {
    stop("scores.transform currently supports only neg_log10", call. = FALSE)
  }
  if (!is.null(config$scores$infinite_score_cap) &&
      (!is.numeric(config$scores$infinite_score_cap) || config$scores$infinite_score_cap <= 0)) {
    stop("scores.infinite_score_cap must be null or a positive number", call. = FALSE)
  }
  assert_positive_integer(config$hhotnet$num_permutations, "hhotnet.num_permutations", allow_zero = TRUE)
  assert_positive_integer(config$hhotnet$lower_size_bound, "hhotnet.lower_size_bound")
  assert_positive_integer(config$hhotnet$cpus, "hhotnet.cpus")
  assert_positive_integer(config$hhotnet$permutation_seed %||% 0L, "hhotnet.permutation_seed", allow_zero = TRUE)
  hhotnet_mode <- config$hhotnet$execution_mode %||% "apptainer"
  if (!hhotnet_mode %in% c("apptainer", "local_python")) {
    stop("hhotnet.execution_mode must be apptainer or local_python", call. = FALSE)
  }
  if (identical(hhotnet_mode, "local_python")) {
    assert_scalar(config$hhotnet$local_dir, "hhotnet.local_dir", "character")
    assert_scalar(config$hhotnet$python_executable %||% "python3", "hhotnet.python_executable", "character")
    if (!is.null(config$hhotnet$compat_dir)) assert_scalar(config$hhotnet$compat_dir, "hhotnet.compat_dir", "character")
  }
  invisible(TRUE)
}

#' Read and validate an analysis1 workflow configuration
#'
#' Relative `project.root` is resolved from the directory containing the YAML.
#' All project input paths are then resolved from that root. Output paths remain
#' under the YAML directory and are guarded by [analysis_paths()].
#'
#' @param file Path to `analysis.yaml`.
#' @param output_root Optional separate output directory for a reproducibility run.
#' @param cache_dir Optional explicit root of cached JEME and HiC resources.
#' @return A validated configuration list with resolved internal metadata.
#' @export
read_analysis_config <- function(file = "analysis.yaml", output_root = NULL, cache_dir = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read analysis.yaml", call. = FALSE)
  }
  file <- normalizePath(file, mustWork = TRUE)
  config <- yaml::read_yaml(file)
  if (isTRUE(config$enhancer_classification$enabled)) {
    config$enhancer_classification$promoter_flank_bp <-
      config$enhancer_classification$promoter_flank_bp %||%
      config$enhancer_classification$gene_node_flank_bp
    config$enhancer_classification$promoter_window <-
      config$enhancer_classification$promoter_window %||% "gene_body"
  }
  validate_analysis_config(config)
  config$config_file <- file
  config$config_dir <- dirname(file)
  config$project_root <- resolve_config_path(config$project$root, config$config_dir)
  config$cache_dir <- resolve_config_path(cache_dir %||% config$resources$cache_dir %||% ".cache",
                                         config$project_root)
  if (!is.null(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    config$config_dir <- normalizePath(output_root, mustWork = TRUE)
    config$project$root <- config$project_root
  }
  config
}

#' Return names of networks available in an analysis
#'
#' `networks.available` defines the variants that can be built. Older
#' configuration files without that field use `networks.selected` as their
#' complete available set.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @return Character vector of unique network names.
#' @export
available_networks <- function(config) {
  if (is.null(config$networks)) stop("config must include networks", call. = FALSE)
  normalise_network_names(config$networks$available %||% config$networks$selected,
                          "networks.available")
}

#' Resolve one or more requested network names
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param networks `NULL` for `networks.selected`, one network name, or a
#'   character vector of names. Every requested name must be in
#'   `networks.available`.
#' @return Character vector in requested order.
#' @export
resolve_networks <- function(config, networks = NULL) {
  available <- available_networks(config)
  requested <- if (is.null(networks)) config$networks$selected else networks
  requested <- normalise_network_names(requested, "networks")
  unknown <- setdiff(requested, available)
  if (length(unknown)) {
    stop("Requested network(s) are not available: ", paste(unknown, collapse = ", "),
         ". Available: ", paste(available, collapse = ", "), call. = FALSE)
  }
  requested
}

#' Return analysis-local output directories
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @return Named character vector of output paths below the YAML directory.
#' @export
analysis_paths <- function(config) {
  required_metadata <- c("config_file", "config_dir", "project_root")
  if (!all(required_metadata %in% names(config))) {
    stop("config must be returned by read_analysis_config()", call. = FALSE)
  }
  root <- normalizePath(config$config_dir, mustWork = TRUE)
  paths <- c(
    root = root,
    logs = file.path(root, "logs"),
    metadata = file.path(root, "metadata"),
    ldak = file.path(root, "LDAK"),
    ldak_annotations = file.path(root, "LDAK", "annotations"),
    ldak_results = file.path(root, "LDAK", "results"),
    ldak_summary = file.path(root, "LDAK", "summary"),
    ldak_scripts = file.path(root, "LDAK", "scripts"),
    hhotnet = file.path(root, "hHotnet"),
    hhotnet_data = file.path(root, "hHotnet", "data"),
    hhotnet_intermediate = file.path(root, "hHotnet", "intermediate"),
    hhotnet_results = file.path(root, "hHotnet", "results"),
    hhotnet_summary = file.path(root, "hHotnet", "summary"),
    hhotnet_cytoscape = file.path(root, "hHotnet", "cytoscape"),
    hhotnet_scripts = file.path(root, "hHotnet", "scripts")
  )
  path_names <- names(paths)
  paths <- stats::setNames(normalizePath(paths, mustWork = FALSE), path_names)
  root_prefix <- paste0(root, .Platform$file.sep)
  outside <- names(paths)[paths != root & !startsWith(paths, root_prefix)]
  if (length(outside)) {
    stop("Refusing output path(s) outside analysis directory: ",
         paste(outside, collapse = ", "), call. = FALSE)
  }
  paths
}

input_paths <- function(config) {
  root <- config$project_root
  paths <- c(
    gwas_summary = resolve_config_path(config$gwas$summary_file, root),
    gwas_pvalues = resolve_config_path(config$gwas$pvalue_file, root),
    gwas_extract = resolve_config_path(config$gwas$extract_file, root),
    ldak_exe = resolve_config_path(config$reference$ldak_exe, root),
    bfile_bed = paste0(resolve_config_path(config$reference$bfile_prefix, root), ".bed"),
    bfile_bim = paste0(resolve_config_path(config$reference$bfile_prefix, root), ".bim"),
    bfile_fam = paste0(resolve_config_path(config$reference$bfile_prefix, root), ".fam")
  )
  if (isTRUE(config$aracne$enabled %||% TRUE)) {
    paths <- c(paths, aracne = resolve_config_path(config$aracne$file, root))
  }
  paths
}

#' Write the fully resolved effective analysis configuration
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param file Optional destination. By default it is written to
#'   `metadata/effective_analysis.yaml` below the analysis directory.
#' @return Normalized path to the written YAML file.
#' @export
write_effective_config <- function(config, file = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to write effective configuration", call. = FALSE)
  }
  paths <- analysis_paths(config)
  if (is.null(file)) file <- file.path(paths[["metadata"]], "effective_analysis.yaml")
  file <- normalizePath(file, mustWork = FALSE)
  root_prefix <- paste0(paths[["root"]], .Platform$file.sep)
  if (!startsWith(file, root_prefix)) {
    stop("Effective configuration must be written below the analysis directory", call. = FALSE)
  }
  effective <- config
  effective$project$root <- config$project_root
  effective$networks$available <- available_networks(config)
  effective$networks$selected <- resolve_networks(config)
  effective$resolved_paths <- list(
    inputs = as.list(input_paths(config)),
    outputs = as.list(paths)
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(effective, file)
  file
}
