#' Internal helper to determine and create a cache directory
#'
#' @param sub_dir Optional sub-directory name to create within the cache.
#' @param cache_dir Optional user-provided base cache directory.
#' @return The absolute path to the (sub-)cache directory.
#' @keywords internal
get_cache_dir <- function(sub_dir = NULL, cache_dir = NULL) {
  cache_dir <- cache_dir %||% getOption("hotnetR2.cache_dir", ".cache")
  if (!is.null(sub_dir)) cache_dir <- file.path(cache_dir, sub_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(cache_dir, mustWork = TRUE)
}

#' Internal helper to determine and create a data_refs directory
#'
#' @param data_refs_dir Optional user-provided directory.
#' @return The absolute path to the data_refs directory.
#' @keywords internal
get_data_refs_dir <- function(data_refs_dir = NULL) {
  if (is.null(data_refs_dir)) {
    # Find project root (reusing logic from get_cache_dir)
    proj_root <- NULL
    if (requireNamespace("here", quietly = TRUE)) {
      proj_root <- tryCatch(here::here(), error = function(e) NULL)
    }

    if (is.null(proj_root)) {
      if (file.exists(".here") || any(grepl("\\.Rproj$", list.files()))) {
        proj_root <- "."
      }
    }

    proposed_refs <- if (!is.null(proj_root)) file.path(proj_root, "data_refs") else "data_refs"

    if (!dir.exists(proposed_refs)) {
      dir.create(proposed_refs, recursive = TRUE)
    }
    data_refs_dir <- proposed_refs
  }
  return(normalizePath(data_refs_dir, mustWork = TRUE))
}

#' Download a file to a cache directory with interactive overwrite prompt
#'
#' @param url The URL to download.
#' @param dest_file The local destination path.
#' @param overwrite Logical, whether to overwrite without prompting.
#' @return The path to the downloaded file.
#' @keywords internal
download_file_cached <- function(url, dest_file, overwrite = FALSE) {
  if (file.exists(dest_file) && !overwrite) {
    if (interactive()) {
      ans <- readline(paste0("File '", dest_file, "' already exists. Overwrite? (y/n): "))
      if (tolower(ans) != "y") {
        message("Using existing file: ", dest_file)
        return(dest_file)
      }
    } else {
      message("File already exists, using: ", dest_file)
      return(dest_file)
    }
  }

  message("Downloading from: ", url)
  utils::download.file(url, destfile = dest_file, mode = "wb")

  if (!file.exists(dest_file) || file.info(dest_file)$size == 0) {
    stop("Failed to download file from ", url)
  }

  return(dest_file)
}

# The following functions are adapted from tmp/r_utils.R

#' Merge preserving the original order if sort==FALSE
#' @keywords internal
merge_with_order <- function(x, y, ..., sort = TRUE) {
  add.id.column.to.data <- function(DATA) {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }

  order.by.id...and.remove.it <- function(DATA) {
    if (!any(colnames(DATA) == "id...")) {
      stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
    }
    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }

  if (!sort) {
    return(order.by.id...and.remove.it(merge(x = add.id.column.to.data(x), y = y, ..., sort = FALSE)))
  } else {
    return(merge(x = x, y = y, ..., sort = sort))
  }
}

#' Download a file when it is not already present
#' @keywords internal
downloadAbsentFile <- function(urlStr, oD = tempdir()) {
  fileName <- utils::tail(strsplit(urlStr, "/")[[1]], 1)
  temp <- file.path(oD, fileName)
  if (!file.exists(temp) || file.info(temp)$size == 0) {
    utils::download.file(urlStr, temp)
  }
  if (file.info(temp)$size == 0) {
    unlink(temp)
    cat(paste0("ERROR: failed to download ", fileName, ".\nPlease check your internet connection."))
    return(NULL)
  }
  return(temp)
}

#' Find rows with NA values in any column
#'
#' @param df A data frame or tibble.
#' @return A data frame containing only rows with at least one NA value.
#' @export
find_na <- function(df) {
  df |> dplyr::filter(dplyr::if_any(dplyr::everything(), is.na))
}

#' Exclude genes located in the MHC region
#'
#' @description
#' This function removes rows from a data frame corresponding to genes
#' located within a specified genomic region (default: the human MHC region
#' on chromosome 6, positions 28,477,797-33,448,354, GRCh38).
#'
#' @param df A data frame containing genomic coordinates. Must include columns:
#'   - `CHR` (chromosome number),
#'   - `START` (start position),
#'   - `END` (end position),
#'   - `gene` (gene identifier).
#' @param chr Chromosome to filter on (default = 6).
#' @param start Start coordinate of the region to exclude (default = 28477797).
#' @param end End coordinate of the region to exclude (default = 33448354).
#' @return A data frame excluding rows within the specified region.
#' @export
exclude_MHC <- function(df, chr = 6, start = 28477797, end = 33448354) {
  idx <- df$CHR == chr & df$START >= start & df$END <= end
  message(paste0("Exclude ", sum(idx), " genes from MHC region"))
  if (sum(idx) > 0) {
    message(paste(df$gene[idx], collapse = ", "))
  }
  df[!idx, ]
}

#' Select rows based on presence of semicolons in a column
#'
#' @param df A data frame or tibble.
#' @param col A column to check for semicolons (unquoted).
#' @param invert Logical (default = FALSE).
#' @return A tibble containing only the rows that satisfy the condition.
#' @export
select_rows_with_semicolon <- function(df, col, invert = FALSE) {
  if (invert) {
    df |> dplyr::filter(!stringr::str_detect({{col}}, ";"))
  } else {
    df |> dplyr::filter(stringr::str_detect({{col}}, ";"))
  }
}

#' Split columns by semicolon into multiple rows
#'
#' @param df A data frame.
#' @param ... One or more unquoted column names (tidy-select) to split.
#' @return A data frame where columns are split into multiple rows.
#' @export
split_columns <- function(df, ...) {
  df |>
    tidyr::separate_rows(..., sep = ";") |>
    dplyr::mutate(dplyr::across(c(...), stringr::str_trim))
}

#' Merge values of a column by group
#'
#' @param df A data frame.
#' @param group_cols Character vector of column names to group by.
#' @param merge_col Name of the column to merge.
#' @param sep Separator to use for merging (default = ";").
#' @return A data frame with merged rows.
#' @export
merge_rows <- function(df, group_cols, merge_col, sep = ";") {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(!!merge_col := paste(unique(.data[[merge_col]]), collapse = sep),
              .groups = "drop")
}

#' Merge two data frames by ENSG IDs (ignoring version suffix)
#'
#' @description This function performs a left join between two data frames based on
#' the main Ensembl gene ID (e.g., ENSG00000141510) only, ignoring any version
#' suffix (e.g., .10). Version numbers in both datasets are automatically
#' stripped before matching to ensure robust harmonization across different annotations.
#'
#' @param df1 A data frame with an `ENSG` column.
#' @param df2 A data frame with an `ENSG` column.
#' @param keep_base Logical, whether to keep only base ENSG IDs (e.g., ENSG00000141510)
#'   in the output. If FALSE, preserves the original versioned ID from `df1`.
#'   Default is TRUE.
#' @return A merged data frame.
#' @export
merge_by_ensg <- function(df1, df2, keep_base = TRUE) {
  df1 <- df1 |> dplyr::mutate(ENSG_base = stringr::str_replace(ENSG, "\\..*", ""))
  df2 <- df2 |> dplyr::mutate(ENSG_base = stringr::str_replace(ENSG, "\\..*", ""))
  merged <- dplyr::left_join(df1, df2, by = "ENSG_base")

  if (keep_base) {
    # Replace original ENSG with version-less base ID
    merged <- merged |>
      dplyr::select(-dplyr::any_of(c("ENSG.x", "ENSG.y"))) |>
      dplyr::rename(ENSG = ENSG_base)
  } else {
    # Preserve the original versioned ID from df1 if it was automatically renamed to ENSG.x
    if ("ENSG.x" %in% colnames(merged)) {
      merged <- merged |>
        dplyr::rename(ENSG = "ENSG.x") |>
        dplyr::select(-dplyr::any_of(c("ENSG.y", "ENSG_base")))
    } else {
      # ENSG should already be there if df2 didn't have an ENSG column before renaming/merging
      merged <- merged |> dplyr::select(-dplyr::any_of(c("ENSG_base")))
    }
  }

  return(merged)
}
