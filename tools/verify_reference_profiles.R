# Regenerate package-owned artifacts in independent directories.
# Usage: Rscript tools/verify_reference_profiles.R PROJECT_ROOT OUTPUT_ROOT
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[[1]])
output <- args[[2]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
profiles <- if (length(args) > 2L) args[-c(1L, 2L)] else c("v2", "v3")
for (profile in profiles) {
  reference <- file.path(project, paste0(profile, "_analysis2"))
  candidate <- file.path(output, profile)
  config <- hotnetR2::read_analysis_config(file.path(reference, "analysis.yaml"),
                                         output_root = candidate)
  hotnetR2::setup_analysis(config, dry_run = FALSE)
  hotnetR2::prepare_ldak_annotations(config, dry_run = FALSE)
  classes <- hotnetR2::classify_enhancers(config)
  print(dplyr::count(classes, .data$enhancer_class))
  files <- c("promoters.loc", "enhancers.loc", "enhancers.details.tsv.gz",
             "enhancer_classification.tsv.gz", "jeme_ensg_harmonization.tsv.gz")
  report <- hotnetR2::compare_artifacts(reference, candidate,
                                       file.path("LDAK", "annotations", files),
                                       decompress = TRUE)
  readr::write_tsv(report, file.path(candidate, "annotation_comparison.tsv"))
  print(report)
  if (!all(report$identical)) stop("Annotation identity failed for ", profile)
  # Raw LDAK results are fixed inputs at this gate. No executable is rerun here.
  for (kind in c("enhancer", "promoter")) {
    subdir <- config$ldak$result_dirs[[kind]]
    file.symlink(file.path(reference, "LDAK", "results", subdir),
                 file.path(candidate, "LDAK", "results", subdir))
  }
  hotnetR2::summarize_analysis_ldak_results(config, dry_run = FALSE)
  scores <- hotnetR2::compare_artifacts(reference, candidate,
    file.path("LDAK/summary", c("enhancer_ldak.tsv.gz", "promoter_ldak.tsv.gz")),
    decompress = TRUE)
  readr::write_tsv(scores, file.path(candidate, "score_comparison.tsv"))
  stopifnot(all(scores$identical))
  hotnetR2::build_hhotnet_networks(config, dry_run = FALSE)
  hotnetR2::validate_hhotnet_inputs(config)
  files <- list.files(file.path(reference, "hHotnet", "data"), "[.]tsv$")
  # Historical source-analysis audits are not current scored network inputs.
  files <- setdiff(files, c("regulatory_edges_class1_pre_harmonization.tsv",
                           "regulatory_edges_class1.tsv"))
  report <- hotnetR2::compare_artifacts(reference, candidate,
                                       file.path("hHotnet", "data", files))
  readr::write_tsv(report, file.path(candidate, "network_comparison.tsv"))
  print(report)
  if (!all(report$identical)) stop("Network identity failed for ", profile)
}
