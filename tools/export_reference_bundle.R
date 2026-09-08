# Export a native-platform archival bundle; never commit the resulting inputs.
# Usage: Rscript tools/export_reference_bundle.R PROJECT_ROOT NEW_BUNDLE_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[[1]])
bundle <- args[[2]]
if (dir.exists(bundle)) stop("Choose a new bundle directory")
config <- hotnetR2::read_analysis_config(file.path(project, "v3_analysis2/analysis.yaml"))
harmonization <- config$regulatory$promoter_harmonization
relative <- c(config$gwas$summary_file, config$gwas$pvalue_file,
  config$gwas$extract_file, paste0(config$reference$bfile_prefix, c(".bed", ".bim", ".fam")),
  config$reference$ncbi_feature_table, harmonization$gtf_file, harmonization$hic_alias_file,
  file.path(".cache/JEME/lasso", list.files(file.path(project, ".cache/JEME/lasso"), "[.]rds$")),
  file.path(".cache/HiC", c("Jung_HiC_P-O.rds", "Jung_HiC_P-P.rds")))
mapping <- tibble::tibble(path = relative, source = file.path(project, relative),
                           category = "frozen_analysis_input")
resolve <- function(path) if (startsWith(path, "/")) path else file.path(project, path)
hotnet <- resolve(config$hhotnet$local_dir)
hotnet_files <- c("LICENSE.txt", "README.md", file.path("src",
  list.files(file.path(hotnet, "src"), recursive = TRUE)))
hotnet_files <- hotnet_files[!grepl("(__pycache__|[.]pyc$)", hotnet_files)]
mapping <- dplyr::bind_rows(mapping,
  tibble::tibble(path = "tools/ldak", source = resolve(config$reference$ldak_exe), category = "native_binary"),
  tibble::tibble(path = file.path("tools/hierarchical-hotnet", hotnet_files),
    source = file.path(hotnet, hotnet_files), category = "hotnet_source_and_native_extension"))
for (i in seq_len(nrow(mapping))) {
  destination <- file.path(bundle, mapping$path[[i]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(mapping$source[[i]], destination, copy.mode = TRUE))
}
manifest <- hotnetR2::capture_resources(bundle, mapping$path) |>
  dplyr::mutate(category = mapping$category, availability = "user_supplied_archive",
    distribution_review = "required_before_redistribution")
original <- purrr::map_chr(mapping$source, ~ digest::digest(file = .x, algo = "sha256"))
stopifnot(identical(original, manifest$sha256))
readr::write_tsv(manifest, file.path(bundle, "resources.tsv"))
yaml::write_yaml(list(platform = R.version$platform,
  hotnet_commit = trimws(processx::run("git", c("-C", hotnet, "rev-parse", "HEAD"))$stdout),
  warning = "Native LDAK and Python extension must match the destination architecture and Python ABI."),
  file.path(bundle, "bundle-environment.yaml"))
message("Exported and verified ", nrow(manifest), " files in ", normalizePath(bundle))
