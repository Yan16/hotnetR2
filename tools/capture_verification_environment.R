# Record the environment used for the full-data development comparisons.
# Usage: Rscript tools/capture_verification_environment.R PROJECT_ROOT OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
config <- hotnetR2::read_analysis_config(file.path(args[[1]], "v3_analysis2/analysis.yaml"))
output <- args[[2]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
python <- config$hhotnet$python_executable
python_info <- processx::run(python, c("-c",
  "import sys, importlib.metadata as m; print(sys.version); print(*sorted(d.metadata['Name']+'=='+d.version for d in m.distributions()), sep=chr(10))"))$stdout
readr::write_lines(python_info, file.path(output, "python-environment.txt"))
readr::write_lines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
ldak <- config$reference$ldak_exe
if (!startsWith(ldak, "/")) ldak <- file.path(config$project_root, ldak)
hotnet <- normalizePath(file.path(config$project_root, config$hhotnet$local_dir), mustWork = FALSE)
if (startsWith(config$hhotnet$local_dir, "/")) hotnet <- config$hhotnet$local_dir
tool_info <- list(
  ldak_sha256 = digest::digest(file = ldak, algo = "sha256"),
  hotnet_commit = trimws(processx::run("git", c("-C", hotnet, "rev-parse", "HEAD"))$stdout),
  hotnet_worktree = processx::run("git", c("-C", hotnet, "status", "--porcelain"))$stdout,
  bedtools_version = trimws(processx::run("bedtools", "--version")$stdout),
  platform = R.version$platform
)
yaml::write_yaml(tool_info, file.path(output, "tools.yaml"))
packages <- as.data.frame(utils::installed.packages(), stringsAsFactors = FALSE)
readr::write_tsv(dplyr::select(packages, "Package", "Version", "Built"),
                 file.path(output, "r-installed-packages.tsv"))
