# Empty tissue selection disables HiC in both annotation and network stages.
# Do not change get_hic(NULL): the public data-reader API still means all tissues.
analysis_hic_tissues <- function(config) {
  tissues <- trimws(as.character(unlist(config$regulatory$hic$tissue_type)))
  tissues[!is.na(tissues) & nzchar(tissues)]
}

analysis_hic_enabled <- function(config) {
  !identical(config$regulatory$hic$enabled, FALSE) && length(analysis_hic_tissues(config)) > 0L
}

empty_analysis_hic <- function() {
  list(
    PO = tibble::tibble(Interacting_fragment = character(), Promoter = character(),
                        CHR = character(), START = integer(), END = integer(), Tissue_type = character()),
    PP = tibble::tibble(Promoter1 = character(), Promoter2 = character(), Tissue_type = character())
  )
}

load_analysis_hic <- function(config, all_tissues = FALSE) {
  if (!analysis_hic_enabled(config)) return(empty_analysis_hic())
  edge_types <- analysis_hic_edge_types(config)
  tissues <- if (all_tissues) NULL else analysis_hic_tissues(config)
  if (setequal(edge_types, c("PO", "PP"))) {
    result <- get_hic(tissue_type = tissues, cache_dir = config$cache_dir)
  } else {
    result <- empty_analysis_hic()
    for (type in edge_types) {
      result[[type]] <- get_hic_by_tissue(tissue_type = tissues, edge_type = type,
                                         cache_dir = config$cache_dir)
    }
  }
  if (isTRUE(config$regulatory$hic$jeme_overlap$enabled) && "PO" %in% edge_types) {
    result <- filter_hic_jeme_overlap(result, config, all_tissues)
  }
  result
}

analysis_hic_edge_types <- function(config) {
  types <- config$regulatory$hic$edge_types
  if (is.null(types)) return(c("PO", "PP"))
  types <- as.character(unlist(types, use.names = FALSE))
  if (!length(types) || anyNA(types) || any(!types %in% c("PO", "PP")) || anyDuplicated(types)) {
    stop("regulatory.hic.edge_types must select PO, PP, or both without duplicates", call. = FALSE)
  }
  types
}
