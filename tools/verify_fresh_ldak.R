# Execute LDAK in a separate result tree using regenerated package annotations.
# Usage: Rscript tools/verify_fresh_ldak.R PROJECT_ROOT ANNOTATION_RUNS OUTPUT_ROOT PROFILE
args <- commandArgs(trailingOnly = TRUE)
profile <- args[[4]]
stopifnot(profile %in% c("v2", "v3"))
reference <- normalizePath(file.path(args[[1]], paste0(profile, "_analysis2")))
candidate <- file.path(args[[3]], profile)
config <- hotnetR2::read_analysis_config(file.path(reference, "analysis.yaml"),
                                        output_root = candidate)
hotnetR2::setup_analysis(config, dry_run = FALSE)
annotations <- file.path(candidate, "LDAK/annotations")
source <- file.path(args[[2]], profile, "LDAK/annotations")
files <- c("enhancers.loc", "promoters.loc", "enhancers.details.tsv.gz", "promoters.details.tsv.gz")
stopifnot(all(file.copy(file.path(source, files), annotations, overwrite = TRUE)))
hotnetR2::run_ldak(config)
hotnetR2::summarize_analysis_ldak_results(config, dry_run = FALSE)
products <- unlist(lapply(config$ldak$result_dirs, function(d) {
  file.path("LDAK/results", d, c("remls.all", "genes.details"))
}))
report <- hotnetR2::compare_artifacts(reference, candidate, products)
readr::write_tsv(report, file.path(candidate, "fresh_ldak_comparison.tsv"))
print(report)
stopifnot(all(report$identical))
scores <- hotnetR2::compare_artifacts(reference, candidate,
  file.path("LDAK/summary", c("enhancer_ldak.tsv.gz", "promoter_ldak.tsv.gz")), decompress = TRUE)
readr::write_tsv(scores, file.path(candidate, "fresh_score_comparison.tsv"))
stopifnot(all(scores$identical))
