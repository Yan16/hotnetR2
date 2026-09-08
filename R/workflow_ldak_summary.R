# Parse and standardise LDAK region-level REML results.

read_ldak_space_table <- function(file, skip = 0L) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required to read LDAK results", call. = FALSE)
  }
  readr::read_table(file, skip = skip, show_col_types = FALSE, progress = FALSE)
}

assert_unique_key <- function(data, key, label) {
  duplicates <- unique(data[[key]][duplicated(data[[key]])])
  if (length(duplicates)) {
    stop(label, " has duplicate ", key, " value(s), e.g. ",
         paste(utils::head(duplicates, 5L), collapse = ", "), call. = FALSE)
  }
}

ldak_summary_specs <- function(config) {
  paths <- analysis_paths(config)
  list(
    enhancer = list(
      result_dir = file.path(paths[["ldak_results"]], ldak_result_subdir(config, "enhancer")),
      detail_file = file.path(paths[["ldak_annotations"]], "enhancers.details.tsv.gz"),
      output_file = file.path(paths[["ldak_summary"]], "enhancer_ldak.tsv.gz")
    ),
    promoter = list(
      result_dir = file.path(paths[["ldak_results"]], ldak_result_subdir(config, "promoter")),
      detail_file = file.path(paths[["ldak_annotations"]], "promoters.details.tsv.gz"),
      output_file = file.path(paths[["ldak_summary"]], "promoter_ldak.tsv.gz")
    )
  )
}

summarize_one_ldak_result <- function(config, node_type, spec) {
  reml_file <- file.path(spec$result_dir, "remls.all")
  if (!file.exists(reml_file)) stop("Missing ", node_type, " LDAK result: ", reml_file, call. = FALSE)
  if (!file.exists(spec$detail_file)) stop("Missing ", node_type, " annotation details: ", spec$detail_file, call. = FALSE)

  result <- read_ldak_space_table(reml_file)
  required <- c("Gene_Name", config$scores$pvalue_column)
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop("Invalid ", node_type, " LDAK schema; missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  assert_unique_key(result, "Gene_Name", paste0(node_type, " remls.all"))
  if (!"SD" %in% names(result)) {
    result$SD <- if ("SE" %in% names(result)) result$SE else NA_real_
  }

  genes_details_file <- file.path(spec$result_dir, "genes.details")
  if (file.exists(genes_details_file)) {
    genes_details <- read_ldak_space_table(genes_details_file, skip = 2L)
    if (all(c("Gene_Name", "Min_Pvalue") %in% names(genes_details))) {
      genes_details <- genes_details[, c("Gene_Name", "Min_Pvalue"), drop = FALSE]
      assert_unique_key(genes_details, "Gene_Name", paste0(node_type, " genes.details"))
      result <- dplyr::left_join(result, genes_details, by = "Gene_Name")
    } else {
      warning("Could not find Gene_Name and Min_Pvalue in ", genes_details_file, call. = FALSE)
    }
  } else {
    warning("Missing genes.details; Min_Pvalue omitted: ", genes_details_file, call. = FALSE)
  }

  details <- readr::read_tsv(spec$detail_file, show_col_types = FALSE)
  if (!"name" %in% names(details)) {
    stop("Annotation details lacks name column: ", spec$detail_file, call. = FALSE)
  }
  details <- dplyr::rename(details, Gene_Name = "name")
  assert_unique_key(details, "Gene_Name", paste0(node_type, " annotation details"))
  unmatched <- sum(!result$Gene_Name %in% details$Gene_Name)
  flank_value <- ldak_flank_bp(config, node_type)
  summary <- dplyr::left_join(result, details, by = "Gene_Name") |>
    dplyr::mutate(
      FDR = stats::p.adjust(.data[[config$scores$pvalue_column]], method = "fdr"),
      gene = .data$Gene_Name,
      node_type = node_type,
      cohort = config$gwas$dataset,
      flank = flank_value
    ) |>
    dplyr::arrange(.data[[config$scores$pvalue_column]], .data$Gene_Name)
  list(data = summary, unmatched_annotations = unmatched, reml_rows = nrow(result))
}

#' Summarize enhancer and promoter LDAK results
#'
#' Standardized tables retain LDAK's `SE` when present and always provide `SD`
#' for compatibility with the established h-HotNet score tables. Joins require
#' unique identifiers and report any LDAK rows without annotation provenance.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param dry_run If `TRUE`, return expected outputs without reading results.
#' @return Named output paths for dry runs, or a data frame summary when run.
#' @export
summarize_analysis_ldak_results <- function(config, dry_run = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("readr", quietly = TRUE)) {
    stop("LDAK summarization requires dplyr and readr", call. = FALSE)
  }
  specs <- ldak_summary_specs(config)
  output_paths <- vapply(specs, `[[`, character(1), "output_file")
  if (isTRUE(dry_run)) return(output_paths)

  paths <- analysis_paths(config)
  dir.create(paths[["ldak_summary"]], recursive = TRUE, showWarnings = FALSE)
  results <- lapply(names(specs), function(node_type) {
    summarize_one_ldak_result(config, node_type, specs[[node_type]])
  })
  names(results) <- names(specs)
  for (node_type in names(results)) {
    readr::write_tsv(results[[node_type]]$data, specs[[node_type]]$output_file, na = "NA")
  }
  report <- data.frame(
    node_type = names(results),
    reml_rows = vapply(results, `[[`, numeric(1), "reml_rows"),
    unmatched_annotations = vapply(results, `[[`, numeric(1), "unmatched_annotations"),
    output_file = unname(output_paths),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(report, file.path(paths[["ldak_summary"]], "summary_manifest.tsv"))
  report
}

ldak_reference_result_file <- function(config, node_type) {
  reference_dir <- config$validation$reference_ldak_dir %||% NULL
  if (is.null(reference_dir)) return(NA_character_)
  reference_dir <- resolve_config_path(reference_dir, config$project_root)
  file.path(
    reference_dir, "results", node_type, "BEEA_Schroder_1k", "remls.all"
  )
}

ldak_reference_summary_file <- function(config, node_type) {
  reference_dir <- config$validation$reference_ldak_dir %||% NULL
  if (is.null(reference_dir)) return(NA_character_)
  reference_dir <- resolve_config_path(reference_dir, config$project_root)
  file_name <- switch(
    node_type,
    enhancer = "Schroder_enhancer_all_JEME_HiC_1k_LDAK.tsv.gz",
    promoter = "Schroder_NCBI_promoter_JEME_HiC_1k_LDAK.tsv.gz",
    stop("Unsupported LDAK node type: ", node_type, call. = FALSE)
  )
  file.path(reference_dir, "summary", file_name)
}

compare_ldak_results <- function(result_file, reference_file, pvalue_column) {
  observed <- read_ldak_space_table(result_file) |>
    dplyr::select("Gene_Name", observed_p = dplyr::all_of(pvalue_column))
  reference <- read_ldak_space_table(reference_file) |>
    dplyr::select("Gene_Name", reference_p = dplyr::all_of(pvalue_column))
  assert_unique_key(observed, "Gene_Name", "observed LDAK result")
  assert_unique_key(reference, "Gene_Name", "reference LDAK result")
  joined <- dplyr::inner_join(observed, reference, by = "Gene_Name")
  delta <- abs(joined$observed_p - joined$reference_p)
  data.frame(
    reference_file = reference_file,
    observed_rows = nrow(observed),
    reference_rows = nrow(reference),
    shared_ids = nrow(joined),
    observed_only_ids = nrow(observed) - nrow(joined),
    reference_only_ids = nrow(reference) - nrow(joined),
    different_lrt_p_perm = sum(delta != 0, na.rm = TRUE),
    max_abs_lrt_p_perm_difference = if (length(delta)) max(delta, na.rm = TRUE) else NA_real_,
    pearson_lrt_p_perm = if (nrow(joined) > 1L) {
      stats::cor(joined$observed_p, joined$reference_p, use = "complete.obs")
    } else NA_real_,
    identical_lrt_p_perm = length(delta) == nrow(observed) && all(delta == 0),
    stringsAsFactors = FALSE
  )
}

compare_ldak_summary_tables <- function(observed_file, reference_file) {
  observed <- readr::read_tsv(observed_file, show_col_types = FALSE, progress = FALSE)
  reference <- readr::read_tsv(reference_file, show_col_types = FALSE, progress = FALSE)
  if (!"Gene_Name" %in% names(observed) || !"Gene_Name" %in% names(reference)) {
    stop("Standardized LDAK summaries must contain Gene_Name", call. = FALSE)
  }
  assert_unique_key(observed, "Gene_Name", "observed standardized LDAK summary")
  assert_unique_key(reference, "Gene_Name", "reference standardized LDAK summary")
  score_columns <- c(
    "Heritability", "SE", "SD", "Null_Likelihood", "Alt_Likelihood", "LRT_Stat",
    "LRT_P_Raw", "LRT_P_Perm", "Min_Pvalue", "FDR"
  )
  common_numeric <- intersect(
    score_columns,
    intersect(names(observed)[vapply(observed, is.numeric, logical(1))],
              names(reference)[vapply(reference, is.numeric, logical(1))])
  )
  joined <- dplyr::inner_join(observed, reference, by = "Gene_Name", suffix = c(".observed", ".reference"))
  comparison <- lapply(common_numeric, function(column) {
    observed_values <- joined[[paste0(column, ".observed")]]
    reference_values <- joined[[paste0(column, ".reference")]]
    comparable <- !is.na(observed_values) & !is.na(reference_values)
    delta <- abs(observed_values[comparable] - reference_values[comparable])
    data.frame(
      metric = column,
      observed_rows = nrow(observed), reference_rows = nrow(reference), shared_ids = nrow(joined),
      observed_only_ids = nrow(observed) - nrow(joined), reference_only_ids = nrow(reference) - nrow(joined),
      comparable_values = sum(comparable), missingness_mismatches = sum(xor(is.na(observed_values), is.na(reference_values))),
      different_values = sum(delta != 0),
      max_abs_difference = if (length(delta)) max(delta) else NA_real_,
      pearson_correlation = if (sum(comparable) > 1L &&
                                stats::sd(observed_values[comparable]) > 0 &&
                                stats::sd(reference_values[comparable]) > 0) {
        stats::cor(observed_values[comparable], reference_values[comparable])
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(comparison)
}

#' Compare standardized LDAK score tables with the configured BEEA reference
#'
#' This read-only comparison does not run LDAK. It compares shared identifiers
#' and every shared numeric score column in the standardized enhancer and
#' promoter tables, then writes a long-format report below `analysis1/LDAK/`.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @return A long-format comparison data frame.
#' @export
compare_ldak_standardized_results <- function(config) {
  specs <- ldak_summary_specs(config)
  rows <- lapply(names(specs), function(node_type) {
    observed_file <- specs[[node_type]]$output_file
    reference_file <- ldak_reference_summary_file(config, node_type)
    if (!file.exists(observed_file)) stop("Missing standardized ", node_type, " summary: ", observed_file, call. = FALSE)
    if (is.na(reference_file) || !file.exists(reference_file)) {
      stop("Missing configured reference ", node_type, " summary: ", reference_file, call. = FALSE)
    }
    dplyr::mutate(
      compare_ldak_summary_tables(observed_file, reference_file),
      node_type = node_type, observed_file = observed_file, reference_file = reference_file,
      .before = 1L
    )
  })
  report <- dplyr::bind_rows(rows)
  output <- file.path(analysis_paths(config)[["ldak_summary"]], "standardized_score_comparison.tsv")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(report, output)
  report
}

#' Validate raw LDAK results and standardized summaries
#'
#' Checks final result files, annotation membership, duplicate identifiers, and
#' optionally compares `LRT_P_Perm` against the configured read-only 1-kb BEEA
#' reference. A compact report is written below `analysis1/LDAK/summary`.
#'
#' @param config A configuration returned by [read_analysis_config()].
#' @param node_types Regional result types to validate.
#' @return A validation data frame.
#' @export
validate_ldak_results <- function(config, node_types = c("enhancer", "promoter")) {
  node_types <- match.arg(node_types, c("enhancer", "promoter"), several.ok = TRUE)
  specs <- ldak_summary_specs(config)
  specs <- specs[node_types]
  pvalue_column <- config$scores$pvalue_column
  rows <- lapply(names(specs), function(node_type) {
    spec <- specs[[node_type]]
    reml_file <- file.path(spec$result_dir, "remls.all")
    details_file <- spec$detail_file
    if (!file.exists(reml_file)) stop("Missing ", node_type, " remls.all: ", reml_file, call. = FALSE)
    if (!file.exists(details_file)) stop("Missing ", node_type, " details: ", details_file, call. = FALSE)
    raw <- read_ldak_space_table(reml_file)
    required <- c("Gene_Name", "SE", pvalue_column)
    missing <- setdiff(required, names(raw))
    if (length(missing)) stop(node_type, " missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    assert_unique_key(raw, "Gene_Name", paste0(node_type, " remls.all"))
    details <- readr::read_tsv(details_file, show_col_types = FALSE)
    assert_unique_key(details, "name", paste0(node_type, " details"))
    annotated <- sum(raw$Gene_Name %in% details$name)
    reference_file <- ldak_reference_result_file(config, node_type)
    comparison <- if (!is.na(reference_file) && file.exists(reference_file)) {
      compare_ldak_results(reml_file, reference_file, pvalue_column)
    } else {
      data.frame(
        reference_file = NA_character_, observed_rows = NA_integer_, reference_rows = NA_integer_,
        shared_ids = NA_integer_, observed_only_ids = NA_integer_, reference_only_ids = NA_integer_,
        different_lrt_p_perm = NA_integer_, max_abs_lrt_p_perm_difference = NA_real_,
        pearson_lrt_p_perm = NA_real_, identical_lrt_p_perm = NA
      )
    }
    cbind(
      data.frame(
        node_type = node_type, result_file = reml_file, result_rows = nrow(raw),
        annotation_rows = nrow(details), annotated_result_rows = annotated,
        duplicate_gene_names = sum(duplicated(raw$Gene_Name)),
        missing_lrt_p_perm = sum(is.na(raw[[pvalue_column]])),
        stringsAsFactors = FALSE
      ),
      comparison
    )
  })
  report <- dplyr::bind_rows(rows)
  out_dir <- analysis_paths(config)[["ldak_summary"]]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(report, file.path(out_dir, "validation_report.tsv"))
  report
}
