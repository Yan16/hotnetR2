#' Run Post-Installation Data Setup
#'
#' @description
#' This function runs the internal setup script to download required
#' datasets (JEME, HiC, NCBI, GTEx, and HGNC) into the package cache. It can
#' optionally create the GENCODE and curated alias references used for promoter
#' harmonization.
#' The datasets are stored in the project's `.cache/` directory by default.
#'
#' @details
#' The function locates and sources the `setup_data.R` script included in the
#' package installation.
#'
#' @param cache_dir Optional. Directory to cache the downloaded files.
#'   If NULL, uses the default project-local cache.
#' @param overwrite Logical. Whether to overwrite existing files without prompting.
#' @param create_harmonization Logical. Whether to create harmonization reference
#'   files after downloading GENCODE. Defaults to FALSE to preserve existing
#'   setup behavior.
#' @param ... Additional arguments passed to specific download functions
#'   (e.g., \code{type = "elasticnet"} for JEME).
#'
#' @return Invisibly returns a list of setup results.
#' @export
#'
#' @examples
#' \dontrun{
#' # Run the full data setup
#' setup_data()
#'
#' # Setup with specific JEME type
#' setup_data(type = "elasticnet")
#'
#' # Also create promoter harmonization references
#' setup_data(create_harmonization = TRUE)
#' }
setup_data <- function(cache_dir = NULL, overwrite = FALSE,
                       create_harmonization = FALSE, ...) {
  script <- system.file("scripts/setup_data.R", package = "hotnetR2")
  if (script == "") {
    stop("Setup script 'inst/scripts/setup_data.R' not found in package installation.")
  }

  message("Sourcing setup script: ", script)

  # Create a local environment to pass parameters to the sourced script
  # Inherit from the current call frame so package functions remain available
  # even when setup_data() is invoked with  rather than library().
  env <- new.env(parent = environment())
  env$cache_dir <- cache_dir
  env$overwrite <- overwrite
  env$create_harmonization <- create_harmonization
  env$args_list <- list(...)

  source(script, local = env)
  return(invisible(env$setup_results))
}

#' @rdname setup_data
#' @export
setup_hotnet_data <- function(cache_dir = NULL, overwrite = FALSE,
                              create_harmonization = FALSE, ...) {
  .Deprecated("setup_data")
  setup_data(
    cache_dir = cache_dir,
    overwrite = overwrite,
    create_harmonization = create_harmonization,
    ...
  )
}
