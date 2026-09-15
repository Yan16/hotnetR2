packaged_alias_tables <- function() {
  data <- new.env(parent = emptyenv())
  utils::data(list = c("alias_link", "alias_link_nodup"), package = "hotnetR2", envir = data)
  list(alias = data$alias, alias_nodup = data$alias_nodup)
}

#' Write the packaged curated alias references to a cache
#'
#' Writes the September 2026 curated gene-alias datasets without requiring
#' GENCODE downloads or importing a GTF. Existing files are retained unless
#' `overwrite = TRUE`. This does not modify analysis configurations or results.
#'
#' @param cache_dir Output cache directory; NULL uses the standard local cache.
#' @param overwrite Replace existing alias RDS files.
#' @return Invisibly, a data frame with artifact, path and status columns.
#' @export
create_alias_references <- function(cache_dir = NULL, overwrite = FALSE) {
  cache <- get_cache_dir(cache_dir = cache_dir)
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  tables <- packaged_alias_tables()
  paths <- file.path(cache, c("alias_link.rds", "alias_link_nodup.rds"))
  write <- overwrite | !file.exists(paths)
  for (i in which(write)) saveRDS(tables[[i]], paths[[i]])
  invisible(tibble::tibble(artifact = names(tables), path = paths,
    status = dplyr::if_else(write, "generated", "skipped_existing")))
}

# New references use raw label -> approved symbol. Unmarked historical caches
# retain their original approved-symbol -> GENCODE-label interpretation.
alias_symbol_pairs <- function(aliases) {
  if (identical(attr(aliases, "mapping_direction"), "raw_to_approved")) {
    dplyr::transmute(aliases, identifier = .data$gene, canonical = .data$Hsym)
  } else {
    dplyr::transmute(aliases, identifier = .data$Hsym, canonical = .data$gene)
  }
}
