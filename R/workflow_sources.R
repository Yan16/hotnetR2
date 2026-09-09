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
  get_hic(tissue_type = if (all_tissues) NULL else analysis_hic_tissues(config),
          cache_dir = config$cache_dir)
}
