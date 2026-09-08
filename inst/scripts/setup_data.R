#!/usr/bin/env Rscript

#' Post-Installation Setup for hotnetR
#'
#' This script downloads the required JEME, HiC, NCBI, and GTEx/GENCODE datasets
#' into the package's local cache (usually the .cache/ directory in your project).
#'
#' You can run this script after installing the package to ensure all data
#' components are ready for the jeme-ldak-workflow and other analyses.

if (!exists("cache_dir")) cache_dir <- NULL
if (!exists("overwrite")) overwrite <- FALSE
if (!exists("create_harmonization")) create_harmonization <- FALSE
if (!exists("args_list")) args_list <- list()

# Helper to merge defaults with args_list and call a function
call_download <- function(func, defaults) {
  # Filter args_list for names that are in the function's signature
  # or just pass everything and let do.call handle it if the function has ...
  # Since most don't have ..., we will merge defaults and only override if in args_list
  final_args <- defaults
  for (n in names(args_list)) {
      final_args[[n]] <- args_list[[n]]
  }
  # Only keep arguments that are valid for the function
  valid_args <- names(formals(func))
  final_args <- final_args[names(final_args) %in% valid_args]

  do.call(func, final_args)
}

message("--- Starting hotnetR Data Setup ---")

setup_results <- list()
total_steps <- if (isTRUE(create_harmonization)) 6L else 5L

# 1. Download JEME interaction data (Enhancers and Promoters)
message("\n[1/", total_steps, "] Downloading JEME data...")
setup_results$jeme <- call_download(download_JEME, list(type = "lasso", cache_dir = cache_dir, overwrite = overwrite))

# 2. Download HiC interaction data
message("\n[2/", total_steps, "] Downloading HiC data...")
setup_results$hic <- call_download(download_HiC, list(cache_dir = cache_dir, overwrite = overwrite))

# 3. Download NCBI feature table
message("\n[3/", total_steps, "] Downloading NCBI feature table...")
setup_results$ncbi <- call_download(download_NBCI_feature_table, list(cache_dir = cache_dir, overwrite = overwrite))

# 4. Download GTEx/GENCODE GTF and metadata
message("\n[4/", total_steps, "] Downloading GTEx and GENCODE data...")
setup_results$gencode <- call_download(download_gtex_gencode, list(cache_dir = cache_dir, overwrite = overwrite))

# 5. Download HGNC complete set
message("\n[5/", total_steps, "] Downloading HGNC complete set...")
setup_results$hgnc <- call_download(download_hgnc_complete, list(cache_dir = cache_dir, overwrite = overwrite))

if (isTRUE(create_harmonization)) {
  message("\n[6/6] Creating promoter harmonization references...")
  setup_results$harmonization <- create_harmonization_references(
    cache_dir = cache_dir,
    overwrite = overwrite
  )
}

message("\n--- Setup Complete! ---")
message("All files are stored in: ", get_cache_dir(cache_dir = cache_dir))
