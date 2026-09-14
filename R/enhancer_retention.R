# Sources whose non-class1 edges are excluded. An explicit YAML list overrides
# the profile default; an empty list retains all classes from both sources.
enhancer_drop_sources <- function(config) {
  settings <- config$enhancer_classification
  configured <- settings$drop_non_class1_edge_sources
  if (!is.null(configured)) return(as.character(unlist(configured, use.names = FALSE)))
  if (identical(settings$promoter_window, "strand_aware_tss")) "HiC_PO" else c("JEME", "HiC_PO")
}

# Annotation rows may represent both sources. Score a region if either source
# retains it, while filtering network edges independently by their source.
retained_enhancer_annotations <- function(classes, config) {
  retained_sources <- setdiff(c("JEME", "HiC_PO"), enhancer_drop_sources(config))
  retained <- classes$enhancer_class == config$enhancer_classification$class1_label
  for (source in retained_sources) {
    retained <- retained | stringr::str_detect(classes$interaction_type,
      paste0("(^|;)", source, "(;|$)"))
  }
  dplyr::filter(classes, retained)
}
