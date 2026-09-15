#' Create promoter harmonization reference files
#'
#' Creates the cleaned GENCODE gene annotation and curated alias files consumed
#' by promoter and HiC harmonization helpers. The alias source bundled with
#' hotnetR2 is curated by the project and is treated as the authoritative default.
#'
#' @param cache_dir Directory containing the `GTEx/` cache. If NULL, uses the
#'   standard project-local cache.
#' @param overwrite Logical; regenerate existing outputs when TRUE.
#' @param alias_file Optional path to a compatible curated alias workbook. When
#'   NULL, uses the packaged September 2026 alias datasets. An explicit workbook
#'   retains the historical workbook interpretation for compatibility.
#' @return Invisibly, a data frame containing each output path and its status.
#' @export
create_harmonization_references <- function(cache_dir = NULL,
                                             overwrite = FALSE,
                                             alias_file = NULL) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "Package 'rtracklayer' is required. Install it with ",
      "BiocManager::install('rtracklayer').",
      call. = FALSE
    )
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required.", call. = FALSE)
  }

  base_cache <- get_cache_dir(cache_dir = cache_dir)
  gtf_file <- file.path(
    base_cache, "GTEx", "gencode.v39.GRCh38.genes.gtf"
  )
  if (!file.exists(gtf_file) || file.info(gtf_file)$size == 0) {
    stop(
      "GENCODE GTF not found or empty: ", gtf_file,
      ". Run download_gtex_gencode(cache_dir = ...) first.",
      call. = FALSE
    )
  }

  if (!is.null(alias_file) && (!nzchar(alias_file) || !file.exists(alias_file) ||
      file.info(alias_file)$size == 0)) {
    stop(
      "Curated alias workbook not found. Supply alias_file or reinstall ",
      "hotnetR2 with its extdata workbook.",
      call. = FALSE
    )
  }

  outputs <- c(
    gencode = file.path(base_cache, "gencode.v39.GRCh38.genes_nodup.rds"),
    alias = file.path(base_cache, "alias_link.rds"),
    alias_nodup = file.path(base_cache, "alias_link_nodup.rds")
  )
  needs_write <- stats::setNames(overwrite | !file.exists(outputs), names(outputs))
  status <- stats::setNames(
    ifelse(needs_write, "generated", "skipped_existing"),
    names(outputs)
  )

  atomic_save_rds <- function(object, path) {
    tmp <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(object, tmp)
    check <- readRDS(tmp)
    if (!is.data.frame(check) || nrow(check) == 0) {
      stop("Generated reference failed validation: ", path, call. = FALSE)
    }
    if (!file.rename(tmp, path)) {
      if (file.exists(path)) unlink(path)
      if (!file.rename(tmp, path)) {
        stop("Could not install generated reference: ", path, call. = FALSE)
      }
    }
  }

  if (needs_write[["gencode"]]) {
    message("Processing GTF annotation from: ", gtf_file)
    gtf <- as.data.frame(rtracklayer::import(gtf_file))
    required_gtf <- c("type", "gene_name", "seqnames")
    missing_gtf <- setdiff(required_gtf, names(gtf))
    if (length(missing_gtf) > 0) {
      stop(
        "GENCODE GTF is missing required columns: ",
        paste(missing_gtf, collapse = ", "), call. = FALSE
      )
    }
    gtf_nodup <- gtf |>
      dplyr::filter(.data$type == "gene") |>
      dplyr::arrange(.data$gene_name, .data$seqnames) |>
      dplyr::group_by(.data$gene_name) |>
      dplyr::distinct(.data$gene_name, .keep_all = TRUE) |>
      dplyr::ungroup()
    if (nrow(gtf_nodup) == 0 || anyDuplicated(gtf_nodup$gene_name)) {
      stop("Cleaned GENCODE reference is empty or contains duplicate gene_name values.", call. = FALSE)
    }
    atomic_save_rds(gtf_nodup, outputs[["gencode"]])
  }

  if (any(needs_write[c("alias", "alias_nodup")])) {
    if (is.null(alias_file)) {
      tables <- packaged_alias_tables()
      for (name in c("alias", "alias_nodup")) {
        if (needs_write[[name]]) atomic_save_rds(tables[[name]], outputs[[name]])
      }
    } else {
    message("Processing curated alias workbook from: ", alias_file)
    alias_link <- as.data.frame(readxl::read_xlsx(alias_file))
    required_alias <- c("gene", "Hsym")
    missing_alias <- setdiff(required_alias, names(alias_link))
    if (length(missing_alias) > 0) {
      stop(
        "Alias workbook is missing required columns: ",
        paste(missing_alias, collapse = ", "), call. = FALSE
      )
    }
    if (nrow(alias_link) == 0) {
      stop("Alias workbook contains no records.", call. = FALSE)
    }
    alias_link_nodup <- alias_link |>
      dplyr::arrange(.data$gene) |>
      dplyr::group_by(.data$gene) |>
      dplyr::distinct(.data$gene, .keep_all = TRUE) |>
      dplyr::ungroup()
    if (anyDuplicated(alias_link_nodup$gene)) {
      stop("Cleaned alias reference contains duplicate gene values.", call. = FALSE)
    }
    if (needs_write[["alias"]]) {
      atomic_save_rds(alias_link, outputs[["alias"]])
    }
    if (needs_write[["alias_nodup"]]) {
      atomic_save_rds(alias_link_nodup, outputs[["alias_nodup"]])
    }
    }
  }

  result <- data.frame(
    artifact = names(outputs),
    path = unname(outputs),
    status = unname(status),
    stringsAsFactors = FALSE
  )
  message("Harmonization reference setup complete in: ", base_cache)
  invisible(result)
}
