#' Initialize a portable v2 or v3 analysis project
#'
#' Creates an analysis YAML using project-relative input paths. Large datasets
#' and external executables are provisioned separately. The profiles preserve
#' the historical v2/v3 decisions, including interval boundary arithmetic.
#'
#' @param path New analysis directory; must not already contain analysis.yaml.
#' @param profile Either `v2` (gene body) or `v3` (TSS).
#' @param tools Named list overriding `ldak_exe`, `hhotnet_dir`, and `python`.
#' @return Path to the created YAML file.
#' @export
initialize_analysis <- function(path, profile = c("v2", "v3"), tools = list()) {
  profile <- match.arg(profile)
  file <- file.path(path, "analysis.yaml")
  if (file.exists(file)) stop("Analysis configuration already exists: ", file, call. = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  config <- yaml::read_yaml(system.file("extdata", paste0(profile, ".yaml"), package = "hotnetR2"))
  config$analysis$name <- basename(normalizePath(path))
  config$project$root <- "."
  config$resources <- list(cache_dir = ".cache", profile = profile)
  config$reference$ldak_exe <- tools$ldak_exe %||% "tools/ldak"
  config$hhotnet$local_dir <- tools$hhotnet_dir %||% "tools/hierarchical-hotnet"
  config$hhotnet$python_executable <- tools$python %||% "python3"
  config$hhotnet$compat_dir <- "tools/compat"
  config$enhancer_classification$source_analysis <- "."
  config$validation <- NULL
  tool_root <- file.path(path, "tools")
  dir.create(tool_root, showWarnings = FALSE)
  file.copy(system.file("compat", package = "hotnetR2"), tool_root, recursive = TRUE)
  yaml::write_yaml(config, file)
  normalizePath(file)
}

#' Capture checksums for a portable input bundle
#'
#' @param root Root containing the input files.
#' @param files Character vector of paths relative to root.
#' @param file Optional TSV manifest destination.
#' @return Tibble with path, byte size and SHA-256 for every input.
#' @export
capture_resources <- function(root, files, file = NULL) {
  validate_resource_paths(files)
  manifest <- tibble::tibble(path = files) |>
    dplyr::mutate(
      size = file.info(file.path(root, .data$path))$size,
      sha256 = purrr::map_chr(.data$path, function(path) {
        digest::digest(file = file.path(root, path), algo = "sha256")
      })
    )
  if (!is.null(file)) readr::write_tsv(manifest, file)
  manifest
}

validate_resource_paths <- function(paths) {
  if (anyNA(paths) || any(!nzchar(paths)) || anyDuplicated(paths) ||
      any(grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$))", paths))) {
    stop("Resource paths must be unique relative paths without parent traversal", call. = FALSE)
  }
}

#' Restore checksummed resources from a local bundle or URL manifest
#'
#' @param manifest Manifest tibble, or path to a TSV with path and sha256 columns.
#'   Remote resources also require a url column.
#' @param root Destination project root.
#' @param source_root Local bundle root, or NULL to use manifest URLs.
#' @param overwrite Whether to replace existing files with different hashes.
#' @return Invisibly, the verified manifest. Files are copied/downloaded into root.
#' @export
provision_resources <- function(manifest, root, source_root = NULL, overwrite = FALSE) {
  if (is.character(manifest)) manifest <- readr::read_tsv(manifest, show_col_types = FALSE)
  if (!all(c("path", "sha256") %in% names(manifest))) stop("Manifest requires path and sha256", call. = FALSE)
  validate_resource_paths(manifest$path)
  if (is.null(source_root) && !"url" %in% names(manifest)) stop("Supply source_root or manifest URLs", call. = FALSE)
  for (i in seq_len(nrow(manifest))) {
    destination <- file.path(root, manifest$path[[i]])
    if (file.exists(destination)) {
      if (identical(digest::digest(file = destination, algo = "sha256"), manifest$sha256[[i]])) next
      if (!overwrite) stop("Existing resource checksum differs: ", destination, call. = FALSE)
    }
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    staged <- tempfile("resource-", tmpdir = dirname(destination))
    if (is.null(source_root)) {
      utils::download.file(manifest$url[[i]], staged, mode = "wb", quiet = TRUE)
    } else {
      if (!file.copy(file.path(source_root, manifest$path[[i]]), staged)) stop("Cannot copy resource: ", manifest$path[[i]], call. = FALSE)
    }
    if (!identical(digest::digest(file = staged, algo = "sha256"), manifest$sha256[[i]])) {
      unlink(staged)
      stop("Resource checksum verification failed: ", manifest$path[[i]], call. = FALSE)
    }
    if (!file.rename(staged, destination)) stop("Cannot install resource: ", destination, call. = FALSE)
  }
  invisible(manifest)
}

#' Compare a declared collection of reproducibility artifacts
#'
#' @param reference Reference analysis directory.
#' @param candidate Independently generated analysis directory.
#' @param files Relative paths of all artifacts required by the comparison.
#' @param decompress Whether gzip files should be compared after decompression.
#'   Defaults to FALSE, requiring raw-byte equality.
#' @return Tibble of paths, reference/candidate hashes and equality indicators.
#' @export
compare_artifacts <- function(reference, candidate, files, decompress = FALSE) {
  hash <- function(root, path) {
    path <- file.path(root, path)
    if (!file.exists(path)) return(NA_character_)
    if (decompress && endsWith(path, ".gz")) {
      return(digest::digest(readr::read_file(path), algo = "sha256", serialize = FALSE))
    }
    digest::digest(file = path, algo = "sha256")
  }
  tibble::tibble(path = files) |>
    dplyr::mutate(reference_sha256 = purrr::map_chr(.data$path, ~ hash(reference, .x)),
                  candidate_sha256 = purrr::map_chr(.data$path, ~ hash(candidate, .x)),
                  identical = !is.na(.data$reference_sha256) & !is.na(.data$candidate_sha256) &
                    .data$reference_sha256 == .data$candidate_sha256)
}
