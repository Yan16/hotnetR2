#!/usr/bin/env Rscript

# Development audit only. Parse source without loading either package.
# Usage: Rscript tools/capture_source_inventory.R PROJECT_ROOT HOTNETR_SOURCE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Supply PROJECT_ROOT and HOTNETR_SOURCE", call. = FALSE)
}
project_root <- normalizePath(args[[1]], mustWork = TRUE)
source_root <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- file.path(project_root, "hotnetR2", "docs", "inventory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sources <- tibble::tibble(
  source = c("hotnetR", "v2_helper", "v3_helper"),
  root = c(source_root,
           file.path(project_root, "v2_analysis2", "hotnethelper"),
           file.path(project_root, "v3_analysis2", "hotnethelper"))
)

function_rows <- function(path, source, relative_path) {
  expressions <- parse(path, keep.source = FALSE)
  purrr::map_dfr(as.list(expressions), function(expr) {
    is_function <- is.call(expr) && is.symbol(expr[[1]]) &&
      as.character(expr[[1]]) %in% c("<-", "=") &&
      length(expr) == 3L && is.symbol(expr[[2]]) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))
    if (!is_function) return(tibble::tibble())
    tibble::tibble(
      source = source,
      file = relative_path,
      name = as.character(expr[[2]]),
      formals = paste(deparse(expr[[3]][[2]], width.cutoff = 500L), collapse = " "),
      body_sha256 = digest::digest(expr[[3]][[3]], algo = "sha256")
    )
  })
}

file_rows <- purrr::pmap_dfr(sources, function(source, root) {
  r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", recursive = TRUE)
  relative_paths <- c("DESCRIPTION", "NAMESPACE", file.path("R", r_files))
  tibble::tibble(source = source, file = relative_paths) |>
    dplyr::mutate(
      sha256 = purrr::map_chr(.data$file, function(path) {
        digest::digest(file = file.path(root, path), algo = "sha256")
      })
    )
})

functions <- purrr::pmap_dfr(sources, function(source, root) {
  files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  exports <- stringr::str_match(namespace, '^export\\(["\']?([^"\')]+)["\']?\\)$')[, 2]
  purrr::map_dfr(files, function(path) {
    function_rows(path, source, file.path("R", basename(path)))
  }) |>
    dplyr::mutate(exported = .data$name %in% exports)
}) |>
  dplyr::arrange(.data$name, .data$source, .data$file)

duplicates <- functions |>
  dplyr::add_count(.data$source, .data$name, name = "definitions_in_source") |>
  dplyr::filter(.data$definitions_in_source > 1L)

collisions <- functions |>
  dplyr::group_by(.data$name) |>
  dplyr::filter(dplyr::n_distinct(.data$source) > 1L,
                dplyr::n_distinct(.data$formals) > 1L) |>
  dplyr::ungroup()

readr::write_tsv(file_rows, file.path(output_dir, "source_files.tsv"))
readr::write_tsv(functions, file.path(output_dir, "function_signatures.tsv"))
readr::write_tsv(duplicates, file.path(output_dir, "duplicate_definitions.tsv"))
readr::write_tsv(collisions, file.path(output_dir, "signature_conflicts.tsv"))
message("Recorded ", nrow(file_rows), " source files and ", nrow(functions),
        " function definitions in ", output_dir)
