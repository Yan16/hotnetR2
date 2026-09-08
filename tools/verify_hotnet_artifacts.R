# Compare freshly executed HotNet artifacts, retaining raw and normalized checks.
# Usage: Rscript tools/verify_hotnet_artifacts.R PROJECT_ROOT CANDIDATE_ROOT
args <- commandArgs(trailingOnly = TRUE)
for (profile in c("v2", "v3")) {
  reference <- normalizePath(file.path(args[[1]], paste0(profile, "_analysis2")))
  candidate <- normalizePath(file.path(args[[2]], profile))
  directories <- c("hHotnet/results", "hHotnet/cytoscape", "hHotnet/summary/clusters")
  inventory <- function(root) unlist(lapply(directories, function(d) {
    file.path(d, list.files(file.path(root, d), recursive = TRUE,
                            pattern = "[.](tsv|graphml|json|cx2)$"))
  }))
  files <- sort(union(inventory(reference), inventory(candidate)))
  report <- hotnetR2::compare_artifacts(reference, candidate, files)
  report$normalization <- dplyr::case_when(
    basename(files) %in% c("export_manifest.tsv", "cluster_table_manifest.tsv") ~ "analysis_root",
    grepl("[.]cyjs[.]json$", files) ~ "producer_name",
    TRUE ~ "none"
  )
  normalize <- function(root, path, rule) {
    file <- file.path(root, path)
    if (!file.exists(file)) return(NA_character_)
    content <- readr::read_file(file)
    if (rule == "analysis_root") content <- stringr::str_replace_all(content,
      stringr::fixed(root), "<ANALYSIS_ROOT>")
    if (rule == "producer_name") content <- stringr::str_replace_all(content,
      stringr::fixed('"generated_by": "hotnethelper::export_hhotnet_graphs"'),
      '"generated_by": "hotnetR2::export_hhotnet_graphs"')
    digest::digest(content, algo = "sha256", serialize = FALSE)
  }
  report$normalized_reference_sha256 <- purrr::map2_chr(files, report$normalization,
    ~ normalize(reference, .x, .y))
  report$normalized_candidate_sha256 <- purrr::map2_chr(files, report$normalization,
    ~ normalize(candidate, .x, .y))
  report$normalized_identical <- !is.na(report$normalized_reference_sha256) &
    !is.na(report$normalized_candidate_sha256) &
    report$normalized_reference_sha256 == report$normalized_candidate_sha256
  readr::write_tsv(report, file.path(candidate, "hotnet_comparison.tsv"))
  print(dplyr::count(report, .data$identical, .data$normalized_identical))
  stopifnot(all(report$normalized_identical))
}
