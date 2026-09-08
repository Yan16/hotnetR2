#' Get HiC data by tissue type
#'
#' @param tissue_type Character vectors of tissue types to retain (e.g., "Esophagus", "Gastric"). If NULL, returns all tissues.
#' @param edge_type Character. Either "PO" (Promoter-Other) or "PP" (Promoter-Promoter)
#' @param cache_dir Root directory containing HiC resources.
#' @return A data frame of the processed HiC interactions filtered to the requested tissue(s)
#' @importFrom readr read_rds
#' @importFrom dplyr filter
get_hic_by_tissue <- function(tissue_type = NULL, edge_type = c("PO", "PP"), cache_dir = NULL) {
  edge_type <- match.arg(edge_type)

  # Check project cache directory (.cache/HiC)
  cache_dir <- get_cache_dir(sub_dir = "HiC", cache_dir = cache_dir)
  file_name <- paste0("Jung_HiC_P-", substr(edge_type, 2, 2), ".rds")
  expected_file <- file.path(cache_dir, file_name)

  if (!file.exists(expected_file)) {
    stop(
      "HiC file '", file_name, "' not found in cache.\n\n",
      "Have you downloaded the HiC dataset yet?\n",
      "Please run the download_HiC function first:\n",
      "  download_HiC()\n",
      call. = FALSE
    )
  }

  dd <- readr::read_rds(expected_file)

  if (!is.null(tissue_type) && length(tissue_type) > 0) {
    tissue_type_lower <- tolower(tissue_type)
    dd <- dplyr::filter(dd, tolower(.data$Tissue_type) %in% tissue_type_lower)
  }

  return(dd)
}
# Example:
# get_hic_by_tissue(c("Esophagus", "Gastric")) |> group_by(Tissue_type) |> tally()

#' Get both PP and PO HiC data by tissue type
#'
#' @param tissue_type Character vectors of tissue types to retain (e.g., "Esophagus", "Gastric"). If NULL, returns all tissues.
#' @return A named list containing 'PP' and 'PO' data frames.
#' @param cache_dir Root directory containing HiC resources.
#' @export
get_hic <- function(tissue_type = NULL, cache_dir = NULL) {
  res <- list(
    PP = get_hic_by_tissue(tissue_type, edge_type = "PP", cache_dir = cache_dir),
    PO = get_hic_by_tissue(tissue_type, edge_type = "PO", cache_dir = cache_dir)
  )
  return(res)
}

# example:
# get_hic("Esophagus")
# get_hic(c("Esophagus", "Gastric"))

# hic <- get_hic()  # for all tissues
# hic$PO |> group_by(Tissue_type) |> tally()
# hic$PP |> group_by(Tissue_type) |> tally()

#' Extract unique enhancer and promoter nodes from HiC data
#'
#' @param hic_data A list containing 'PP' and 'PO' data frames returned by get_hic()
#' @return A named list with two elements: `enhancers` and `promoters`, each a data frame of unique regions with an aggregated `tissues` column
#' @importFrom dplyr select all_of distinct group_by across summarise bind_rows .data
#' @export
extract_hic_nodes <- function(hic_data) {
  # PO: Interacting_fragment (Enhancer), Promoter (Target)
  # PP: Promoter1, Promoter2 (both Promoters)

  # 1. Process PO for Enhancers and Promoters
  po <- hic_data$PO
  if (!is.null(po) && nrow(po) > 0) {
    # Enhancers from PO
    # Interacting_fragment is already parsed into CHR, START, END in download_HiC
    # but let's make sure columns exist
    enh_cols <- c("Interacting_fragment", "CHR", "START", "END")
    available_enh <- enh_cols[enh_cols %in% colnames(po)]

    enhancers <- po |>
      dplyr::select(dplyr::all_of(c(available_enh, "Tissue_type"))) |>
      dplyr::distinct() |>
      dplyr::group_by(dplyr::across(dplyr::all_of(available_enh))) |>
      dplyr::summarise(
        tissues = paste(unique(.data$Tissue_type), collapse = ", "),
        .groups = "drop"
      )

    # Extract Promoters from PO
    # Need to split if multi-gene
    promoters_po <- po |>
      split_columns(Promoter) |>
      dplyr::select(Promoter, Tissue_type) |>
      dplyr::distinct()
  } else {
    enhancers <- data.frame()
    promoters_po <- data.frame()
  }

  # 2. Process PP for Promoters
  pp <- hic_data$PP
  if (!is.null(pp) && nrow(pp) > 0) {
    promoters_pp1 <- pp |>
      split_columns(Promoter1) |>
      dplyr::select(Promoter = Promoter1, Tissue_type) |>
      dplyr::distinct()

    promoters_pp2 <- pp |>
      split_columns(Promoter2) |>
      dplyr::select(Promoter = Promoter2, Tissue_type) |>
      dplyr::distinct()

    promoters_pp <- dplyr::bind_rows(promoters_pp1, promoters_pp2) |>
      dplyr::distinct()
  } else {
    promoters_pp <- data.frame()
  }

  # 3. Combine and aggregate Promoters
  all_promoters <- dplyr::bind_rows(promoters_po, promoters_pp)

  if (nrow(all_promoters) > 0) {
    promoters <- all_promoters |>
      dplyr::group_by(Promoter) |>
      dplyr::summarise(
        tissues = paste(unique(.data$Tissue_type), collapse = ", "),
        .groups = "drop"
      )
  } else {
    promoters <- data.frame()
  }

  # 4. Final NA filtering
  if (nrow(enhancers) > 0) {
    enhancers <- enhancers |>
      dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
  }
  if (nrow(promoters) > 0) {
    promoters <- promoters |>
      dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
  }

  list(enhancers = enhancers, promoters = promoters)
}

#' Harmonize HiC networks using preprocessed alias link
#'
#' @param hic A list containing PO and PP networks
#' @param alias_link_file Path to alias link RDS
#' @param verbose Logical.
#' @param rm_na Logical.
#' @return Harmonized HiC networks
#' @export
harmonize_hic_gtf <- function(hic, alias_link_file, verbose = TRUE, rm_na = FALSE) {
  # Load alias link
  if (verbose) message("Load alias link file: ", alias_link_file)
  alias_link <- readr::read_rds(alias_link_file) |>
    dplyr::select(gene, Hsym) |>
    dplyr::filter(!is.na(Hsym)) |>
    dplyr::ungroup()

  # Harmonize PO
  this_PO <- hic$PO |>
    dplyr::left_join(alias_link |> dplyr::rename(Target_harmonized = gene, Target = Hsym), by = "Target") |>
    dplyr::mutate(Target_harmonized = dplyr::if_else(is.na(Target_harmonized), Target, Target_harmonized))

  if (nrow(find_na(this_PO)) > 0 && rm_na) {
    this_PO <- this_PO |> dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
  }

  # Harmonize PP
  this_PP <- hic$PP |>
    dplyr::left_join(alias_link |> dplyr::rename(Promoter1_harm = gene, Promoter1 = Hsym), by = "Promoter1") |>
    dplyr::mutate(Promoter1_harm = dplyr::if_else(is.na(Promoter1_harm), Promoter1, Promoter1_harm)) |>
    dplyr::left_join(alias_link |> dplyr::rename(Promoter2_harm = gene, Promoter2 = Hsym), by = "Promoter2") |>
    dplyr::mutate(Promoter2_harm = dplyr::if_else(is.na(Promoter2_harm), Promoter2, Promoter2_harm))

  if (nrow(find_na(this_PP)) > 0 && rm_na) {
    this_PP <- this_PP |> dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
  }

  return(list(PO = this_PO, PP = this_PP))
}

#' Deprecated: Load HiC networks
#'
#' @description This function is deprecated. Please use \code{\link{get_hic}} instead.
#' @param param Parameters list.
#' @param verbose Logical.
#' @param rm_na Logical.
load_hic <- function(param, verbose = TRUE, rm_na = FALSE) {
  .Deprecated("get_hic", package = "hotnetR2",
             msg = "load_hic() is legacy. Please switch to get_hic() for modern data access.")
  get_hic(tissue_type = param$HiC_param$tissue)
}
