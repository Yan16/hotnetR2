#' Get JEME data by tissue nfile
#'
#' @param nfile Numeric or character. The tissue numeric file identifier (e.g., 77, 92)
#' @param method Character. Either "lasso" or "elasticnet"
#' @param cache_dir Root directory containing JEME and HiC resources.
#' @return A data frame of the processed JEME interactions for that tissue
#' @importFrom readr read_rds
get_jeme_by_tissue <- function(nfile, method = c("lasso", "elasticnet"), cache_dir = NULL) {
  method <- match.arg(method)
  nfile <- as.character(nfile)

  # Check project cache directory (.cache/JEME)
  cache_dir <- get_cache_dir(sub_dir = "JEME", cache_dir = cache_dir)
  sys_type_dir <- file.path(cache_dir, method)
  expected_file <- file.path(sys_type_dir, paste0("encoderoadmap_", method, ".", nfile, ".rds"))

  if (!file.exists(expected_file)) {
    stop(
      "JEME file 'encoderoadmap_", method, ".", nfile, ".rds' not found in cache.\n\n",
      "Have you downloaded the JEME dataset yet?\n",
      "Please run the download_JEME function first:\n",
      "  download_JEME(type = \"", method, "\")\n",
      call. = FALSE
    )
  }

  readr::read_rds(expected_file) |> tibble::add_column(nfile = nfile)
}
# example
# jeme_77 <- get_jeme_by_tissue(77, method = "lasso")


#' Get JEME data for specific tissues
#'
#' @param key_column Character. The column name in the key_tissues_jeme dataset to match (e.g., "desc1"). If NULL (default), reads all tissues.
#' @param value Character vector. The values to match in the key_column.
#' @param method Character. Either "lasso" or "elasticnet"
#' @param simplified Logical. If TRUE, concatenates the resulting list into a single data frame
#' @param cache_dir Root directory containing the JEME resource directory.
#' @return A list of data frames, or a single combined data frame if simplified=TRUE
#' @importFrom purrr map
#' @importFrom dplyr filter bind_rows left_join select distinct group_by summarise across all_of any_of
#' @importFrom stats setNames
#' @export
get_jeme <- function(key_column = NULL, value = NULL, method = c("lasso", "elasticnet"), simplified = TRUE, cache_dir = NULL) {
  method <- match.arg(method)

  # Load the internal key dataset
  env <- new.env()
  utils::data("key_tissues_jeme", package = "hotnetR2", envir = env)
  keys <- env$key_tissues_jeme

  if (is.null(keys)) {
    # Fallback to loading it if we are developing locally
    if (file.exists("data/key_tissues_jeme.rda")) {
      load("data/key_tissues_jeme.rda", envir = env)
      keys <- env$key_tissues_jeme
    } else {
      stop("Internal dataset key_tissues_jeme not found.")
    }
  }

  # Ensure we only process files for the requested method to avoid lasso/elasticnet duplication
  if ("file" %in% colnames(keys)) {
    keys <- keys[grepl(method, keys$file, ignore.case = TRUE), , drop = FALSE]
  }

  if (!is.null(key_column) && !is.null(value)) {
    if (!key_column %in% colnames(keys)) {
      stop("Column '", key_column, "' not found in key_tissues_jeme dataset.")
    }

    # Check for missing values case-insensitively
    key_col_vals <- tolower(as.character(keys[[key_column]]))
    val_check <- tolower(as.character(value))
    missing_vals <- value[!val_check %in% key_col_vals]

    if (length(missing_vals) > 0) {
      warning(
        "The following values were not found in column '", key_column, "': ",
        paste(missing_vals, collapse = ", "), ". They may be misspelled."
      )
    }

    keys <- keys[key_col_vals %in% val_check, , drop = FALSE]
  }

  if (nrow(keys) == 0) {
    stop("No tissues found matching the specified criteria.")
  }

  list_names <- keys$nfile
  if ("file" %in% colnames(keys)) {
    list_names <- keys$file
  }

  res_list <- keys$nfile |>
    stats::setNames(list_names) |>
    purrr::map(~ get_jeme_by_tissue(.x, method = method, cache_dir = cache_dir))

  if (simplified) {
    cmb <- dplyr::bind_rows(res_list, .id = "file")
    if ("desc1" %in% colnames(keys) && "file" %in% colnames(keys)) {
      cmb <- dplyr::left_join(cmb, keys[, c("file", "desc1")], by = "file")
    }
    return(cmb)
  }
  return(res_list)
}
# ex:
# jeme_data <- get_jeme(key_column = "nfile", value = c("77", "90", "92", "108"), method = "lasso")

# ex: shpuld be same as above
# jeme_data <- get_jeme(key_column = "cat2", value = c("GI_ESOPHAGUS", GI_STOMACH"), method = "lasso")

#' Extract unique enhancer and promoter nodes from JEME data
#'
#' @param jeme_data A data frame (simplified) or a list of data frames (non-simplified) returned by get_jeme()
#' @return A named list with two elements:
#' \describe{
#'   \item{enhancers}{A data frame of unique enhancer regions with the columns:
#'     \itemize{
#'       \item{`enhancer`}
#'       \item{`CHR`}
#'       \item{`START`}
#'       \item{`END`}
#'       \item{`files`} (aggregated file identifiers where the enhancer appears)
#'       \item{`tissues`} (aggregated tissue descriptions)
#'     }
#'   }
#'   \item{promoters}{A data frame of unique promoter nodes with the columns:
#'     \itemize{
#'       \item{`promoterFull`}
#'       \item{`promoter`}
#'       \item{`CHR2`}
#'       \item{`location`}
#'       \item{`files`} (aggregated file identifiers where the promoter appears)
#'       \item{`tissues`} (aggregated tissue descriptions)
#'     }
#'   }
#' }
#' Each returned data frame contains aggregated identifier columns (`files`, `nfiles`)
#' and a `tissues` column when available. Column names reflect the conventions used
#' by the input JEME data and by `promoter_range()` processing.
#' @importFrom dplyr mutate select arrange distinct case_when across .data all_of group_by summarise filter if_all bind_rows
#' @export
extract_jeme_nodes <- function(jeme_data) {
  if (is.list(jeme_data) && !is.data.frame(jeme_data)) {
    jeme_data <- dplyr::bind_rows(jeme_data, .id = "file")
  }
  # the column names in jeme_data should be
  # [1] "file"         "enhancer"     "CHR"          "START"        "END"
  # [6] "promoterFull" "ENSG"         "promoter"     "CHR2"         "location"
  # [11] "conf_score"   "nfile"        "desc1"
  id_cols <- c("file", "nfile")
  id_cols <- id_cols[id_cols %in% colnames(jeme_data)]

  desc_col <- if ("desc1" %in% colnames(jeme_data)) "desc1" else NULL

  # Internal helper for efficient aggregation
  aggregate_nodes <- function(data, cols, ids, desc) {
    # 1. Select only needed columns and deduplicate rows early for efficiency
    relevant_cols <- c(cols, ids)
    if (!is.null(desc)) relevant_cols <- c(relevant_cols, desc)

    data |>
      dplyr::select(dplyr::all_of(relevant_cols)) |>
      dplyr::distinct() |> # Efficient early deduplication
      dplyr::group_by(dplyr::across(dplyr::all_of(cols))) |>
      dplyr::summarise(
        # Aggregate multiple IDs into pluralized columns (e.g., file -> files, nfile -> nfiles)
        dplyr::across(dplyr::all_of(ids), ~ paste(unique(.x), collapse = ", "), .names = "{.col}s"),
        # Aggregate tissues if desc exists
        tissues = if (!is.null(desc)) paste(unique(.data[[desc]]), collapse = ", ") else NA_character_,
        .groups = "drop"
      ) |>
      # Final filtering for NA nodes (ensures unique keys don't contain NAs)
      dplyr::filter(dplyr::if_all(dplyr::all_of(cols), ~ !is.na(.)))
  }

  # Define core columns for enhancers and promoters
  enhancer_cols <- c("enhancer", "CHR", "START", "END")
  # For promoters, we want standardized coordinate names after ranging
  target_promoter_cols <- c("promoterFull", "promoter", "ENSG", "CHR", "START", "END")

  # Filter to columns available in the input data for enhancers
  available_enhancer <- enhancer_cols[enhancer_cols %in% colnames(jeme_data)]

  # Aggregate if any columns are found
  enhancers <- if (length(available_enhancer) > 0) {
    aggregate_nodes(jeme_data, available_enhancer, id_cols, desc_col)
  } else {
    data.frame()
  }
  enhancers <- enhancers |>
    dplyr::mutate(CHR = sub("^chr", "", as.character(.data$CHR), ignore.case = TRUE),
        CHR = dplyr::case_when(
          .data$CHR == "X" ~ "23",
          .data$CHR == "Y" ~ "24",
          TRUE ~ .data$CHR
        ),
        START = as.integer(.data$START),
        END = as.integer(.data$END)
      ) |> dplyr::arrange(.data$CHR, .data$START, .data$END)

  promoters <- if ("location" %in% colnames(jeme_data)) {
    # 1. Standardize promoter coordinates using the utility
    # We select promoter-related columns to avoid column collisions with enhancer CHR/START/END
    promoter_base_cols <- c("promoterFull", "promoter", "ENSG", "CHR2", "location")
    existing_base <- promoter_base_cols[promoter_base_cols %in% colnames(jeme_data)]

    processed_proms <- jeme_data |>
      dplyr::select(dplyr::all_of(c(existing_base, id_cols, desc_col))) |>
      promoter_range(source = "jeme")

    # 2. Identify available columns for aggregation
    # We include legacy names (CHR2, location) to satisfy downstream expectations
    target_promoter_cols <- c("promoterFull", "promoter", "ENSG", "CHR2", "location", "CHR", "START", "END")
    available_target <- target_promoter_cols[target_promoter_cols %in% colnames(processed_proms)]

    aggregate_nodes(processed_proms, available_target, id_cols, desc_col)
  } else {
    data.frame()
  }

  promoters <- promoters |>
    dplyr::mutate(CHR = sub("^chr", "", as.character(.data$CHR), ignore.case = TRUE),
        CHR = dplyr::case_when(
          .data$CHR == "X" ~ "23",
          .data$CHR == "Y" ~ "24",
          TRUE ~ .data$CHR
        ) |> as.integer(),
        START = as.integer(.data$START),
        END = as.integer(.data$END)
      ) |> dplyr::arrange(.data$CHR, .data$START, .data$END)


  # Clean up 'tissues' column if it's entirely NA (meaning desc_col was NULL)
  if (!is.null(enhancers$tissues) && all(is.na(enhancers$tissues))) enhancers$tissues <- NULL
  if (!is.null(promoters$tissues) && all(is.na(promoters$tissues))) promoters$tissues <- NULL

  list(enhancers = enhancers, promoters = promoters)
}
# ex:
# jeme_nodes <- extract_jeme_nodes(jeme_data)
# lapply(jeme_nodes, head)

#' Define promoter range based on source type
#'
#' This utility converts a 1-based biological TSS to its one-base, 0-based
#' half-open BED interval `[TSS-1,TSS)`.
#'
#' @param x A data frame containing coordinate information. Either NCBI feature table format (with `start` column) or JEME format (with `location` column).
#' @param source Character. Either "NCBI" (expects `start` column) or "jeme" (expects `location`).
#' @return A data frame with normalized CHR, START, and END columns.
#' @importFrom dplyr mutate rename any_of .data
#' @export
promoter_range <- function(x, source = c("NCBI", "jeme")) {
  source <- match.arg(source)

  if (source == "NCBI") {
    if (!"start" %in% colnames(x)) stop("Source df must contain 'start' column for NCBI source.")
    if ("strand" %in% colnames(x) && !"end" %in% colnames(x)) {
      stop("NCBI source with 'strand' must also contain 'end'.")
    }
    tss <- if (all(c("strand", "end") %in% colnames(x))) {
      if (any(!x$strand %in% c("+", "-"))) stop("NCBI source contains an unsupported strand.")
      ifelse(x$strand == "+", as.integer(x$start), as.integer(x$end))
    } else {
      as.integer(x$start)
    }
    res <- x |>
      dplyr::mutate(
        START = tss - 1L,
        END = tss
      )
    if ("chromosome" %in% colnames(res)) {
      res <- dplyr::rename(res, CHR = "chromosome")
    }
  } else {
    if (!"location" %in% colnames(x)) stop("Source df must contain 'location' column for jeme source.")
    res <- x |>
      dplyr::mutate(
        START = as.integer(.data$location) - 1L,
        END = as.integer(.data$location)
      )
    if ("CHR2" %in% colnames(res)) {
      res <- dplyr::rename(res, CHR = "CHR2")
    }
  }

  # Ensure CHR is present if not already renamed (e.g. if already called CHR)
  if (!"CHR" %in% colnames(res) && "chromosome" %in% colnames(res)) {
    res <- dplyr::rename(res, CHR = "chromosome")
  }

  return(res)
}
# ex:
# For NCBI feature table format:
# r1 <- promoter_range(ncbi_df, source = "NCBI") # TODO fix this with real NCBI data
# r2 <- promoter_range(nodes$promoters, source = "jeme")


#' Harmonize gene symbols by Ensembl ID with GENCODE annotation, to match those in GTEx V10
#'
#' @param df A data frame containing at least `ENSG` and `promoter` columns.
#' @param gtf A tibble/dataframe containing GENCODE annotation. If NULL (default),
#'   attempts to read the cached GTEx GENCODE GTF.
#' @param verbose Logical, whether to print progress messages (default = TRUE)
#' @return A data frame with an added `promoter_harmonized` column.
#' @export
harmonize_promoter_by_ensg <- function(df, gtf = NULL, verbose = TRUE) {
  if (verbose) message("Harmonizing promoters by Ensembl ID with GENCODE GTF...")
  if (verbose) message("We match based on the main Ensembl ID (e.g., ENSG00000141510) only, without the gene version (e.g., .10)...")

  if (!all(c("ENSG", "promoter") %in% colnames(df))) {
    stop("Input data frame must contain 'ENSG' and 'promoter' columns.")
  }

  if (is.null(gtf)) {
    if (verbose) message("No GTF provided. Reading GTEx V10 GTF from default cache...")
    gtf <- read_gtex_gencode_gtf()
  }

  # Ensure GTF has expected 'ENSG' and 'gene' columns for merge_by_ensg
  # rtracklayer usually produces 'gene_id' and 'gene_name'
  if (!"ENSG" %in% colnames(gtf) && "gene_id" %in% colnames(gtf)) {
    gtf <- dplyr::rename(gtf, ENSG = "gene_id")
  }
  if (!"gene" %in% colnames(gtf) && "gene_name" %in% colnames(gtf)) {
    gtf <- dplyr::rename(gtf, gene = "gene_name")
  }

  # Preserve original order and use merge_by_ensg
  # We use keep_base = FALSE to ensure original versioned IDs from 'df' are preserved in the output.
  # gtf is renamed to promoter_harmonized to match the expected return
  res_df <- df |>
    dplyr::mutate(original_order... = seq_len(nrow(df))) |>
    merge_by_ensg(gtf |> dplyr::rename(promoter_harmonized = "gene"), keep_base = FALSE) |>
    dplyr::arrange(.data$original_order...)

  # Fallback to original promoter if not found in GTF (NA)
  res_df <- res_df |>
    dplyr::mutate(promoter_harmonized = dplyr::if_else(is.na(.data$promoter_harmonized), .data$promoter, .data$promoter_harmonized)) |>
    dplyr::select(-dplyr::any_of("original_order..."))

  return(res_df)
}

#' Harmonize JEME promoter gene names with GENCODE annotation
#'
#' @description This function is deprecated. Please use \code{\link{harmonize_promoter_by_ensg}} instead.
#' @param jeme A tibble containing JEME enhancer-promoter edges. Must contain `ENSG`.
#' @param gtf A tibble containing GENCODE annotation. If NULL (default),
#'   attempts to read the cached GTEx GENCODE GTF.
#' @param verbose Logical, whether to print progress messages (default = TRUE)
#' @return A tibble with harmonized JEME data
#' @export
harmonize_jeme_genecode <- function(jeme, gtf = NULL, verbose = TRUE) {
  .Deprecated("harmonize_promoter_by_ensg",
    package = "hotnetR2",
    msg = "harmonize_jeme_genecode() is legacy. Please use the more general harmonize_promoter_by_ensg()."
  )

  if (verbose) message("Harmonizing JEME promoters with GENCODE GTF...")

  # Handle legacy 'Target' vs 'promoter' naming
  if (!"promoter" %in% colnames(jeme) && "Target" %in% colnames(jeme)) {
    jeme <- dplyr::rename(jeme, promoter = "Target")
  }

  # Use the new general function
  res <- harmonize_promoter_by_ensg(jeme, gtf = gtf, verbose = verbose)

  # Maintain legacy return column name alongside new one if needed
  # But the user specifically asked for promoter_harmonized
  return(res)
}

#' Deprecated: Summarize harmonized JEME results
#'
#' @description This function is deprecated.
#' @param jeme_harmonized A tibble returned by `harmonize_jeme_genecode()` or `harmonize_promoter_by_ensg()`
#' @param max_print Maximum items to print.
#' @param verbose Logical.
#' @return A summary list.
summary_harmonized_jeme <- function(jeme_harmonized, max_print = 10, verbose = TRUE) {
  .Deprecated(msg = "summary_harmonized_jeme() is legacy and may not work with new harmonization outputs.")

  if (verbose) message("Summarizing harmonized JEME data...")

  # Handle legacy 'Target' vs 'promoter' naming if needed for summary
  target_col <- if ("Target" %in% colnames(jeme_harmonized)) "Target" else "promoter"
  harmonized_col <- if ("Target_harmonized" %in% colnames(jeme_harmonized)) "Target_harmonized" else "promoter_harmonized"

  if (!all(c(target_col, harmonized_col) %in% colnames(jeme_harmonized))) {
    warning("Input data frame missing expected columns for summary.")
    return(NULL)
  }

  summary_tbl <- table(jeme_harmonized[[target_col]] == jeme_harmonized[[harmonized_col]]) |> as.data.frame()
  changed_promoters <- unique(jeme_harmonized[[target_col]][jeme_harmonized[[target_col]] != jeme_harmonized[[harmonized_col]]])

  if (verbose) {
    message("How many promoters were changed (FALSE = changed):")
    print(summary_tbl)
    message("Number of unique changed promoters: ", length(changed_promoters))
    print(utils::head(changed_promoters, max_print))
  }
  return(list(summary_tbl = summary_tbl, changed_promoters = changed_promoters))
}

#' Deprecated: Load JEME network
#'
#' @description This function is deprecated. Please use \code{\link{get_jeme}} instead.
#' @param param Parameters list.
#' @param verbose Logical.
load_jeme <- function(param, verbose = TRUE) {
  .Deprecated("get_jeme",
    package = "hotnetR2",
    msg = "load_jeme() is legacy. Please switch to get_jeme() for modern data access."
  )
  # Wrapper attempt to map old param structure to new call
  get_jeme(key_column = "desc1", value = param$JEME_param$tissue, method = param$JEME_param$method)
}

#' Read and process a single JEME CSV file
#'
#' @param infile Path to the JEME CSV file.
#' @param out_rds Optional path to save the processed RDS file. If NULL (default), does not save.
#' @param simplify_col Logical, whether to drop the `conf_score` column (default FALSE).
#' @return A data frame of processed JEME interactions.
#' @importFrom readr read_csv write_rds
#' @importFrom tidyr separate
#' @importFrom dplyr rename select any_of
#' @export
read_jeme_csv <- function(infile, out_rds = NULL, simplify_col = FALSE) {
  dd <- readr::read_csv(infile, col_names = FALSE, show_col_types = FALSE) |>
    tidyr::separate("X1", c("CHR", "START", "END"), remove = FALSE, convert = TRUE) |>
    tidyr::separate("X2", c("ENSG", "promoter", "CHR2", "location", "strand"), sep = "\\$", remove = FALSE, convert = TRUE) |>
    dplyr::rename(enhancer = "X1", promoterFull = "X2")

  if ("X3" %in% colnames(dd)) {
    dd <- dd |> dplyr::rename(conf_score = "X3")
  }

  # dd <- dd |> dplyr::select(-dplyr::any_of("var1"))

  if (!is.null(out_rds)) {
    if (!dir.exists(dirname(out_rds))) dir.create(dirname(out_rds), recursive = TRUE)
    readr::write_rds(dd, out_rds, compress = "gz")
  }

  if (simplify_col) {
    dd <- dd |> dplyr::select(-dplyr::any_of(c("conf_score")))
  }
  return(dd)
}
