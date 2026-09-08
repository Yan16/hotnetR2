#' Download GTEx TPM data for a specific tissue
#'
#' @param tissue Tissue name (e.g., "Esophagus_Gastroesophageal_Junction", "Whole_Blood")
#' @param datasetId GTEx dataset ID (default "gtex_v10")
#' @param cache_dir Directory to cache the downloaded file. If NULL, proposes .cache/ in R project or current dir.
#' @param overwrite Logical, whether to overwrite existing file without prompting (default FALSE)
#' @return Path to the downloaded file
#' @importFrom gtexr get_tissue_site_detail
#' @importFrom utils download.file
#' @export
download_gtex_tpm <- function(tissue = "Esophagus_Gastroesophageal_Junction", datasetId = "gtex_v10", cache_dir = NULL, overwrite = FALSE) {

  # 1. Determine/Create cache directory (using .cache/GTEx)
  base_cache <- get_cache_dir(cache_dir = cache_dir)
  gtex_cache <- file.path(base_cache, "GTEx")
  if (!dir.exists(gtex_cache)) dir.create(gtex_cache, recursive = TRUE)

  # 2. Prepare URL and local file path
  tissue_lower <- tolower(tissue)
  version <- gsub("gtex_", "", datasetId) # e.g., "v10"
  file_name <- paste0("gene_tpm_", version, "_", tissue_lower, ".gct.gz")
  url <- paste0(
    "https://storage.googleapis.com/adult-gtex/bulk-gex/",
    version,
    "/rna-seq/tpms-by-tissue/", file_name
  )

  dest_file <- file.path(gtex_cache, file_name)

  # 3. Download and return
  return(download_file_cached(url, dest_file, overwrite = overwrite))
}

#' Download GENCODE annotation for GTEx
#'
#' @param cache_dir Directory to cache the downloaded file. If NULL, proposes .cache/ in R project or current dir.
#' @param overwrite Logical, whether to overwrite existing file without prompting (default FALSE)
#' @return Path to the downloaded file
#' @export
download_gtex_gencode <- function(cache_dir = NULL, overwrite = FALSE) {
  base_cache <- get_cache_dir(cache_dir = cache_dir)
  gtex_cache <- file.path(base_cache, "GTEx")
  if (!dir.exists(gtex_cache)) dir.create(gtex_cache, recursive = TRUE)

  url <- "https://storage.googleapis.com/adult-gtex/references/v10/reference-tables/gencode.v39.GRCh38.genes.gtf"
  dest_file <- file.path(gtex_cache, basename(url))
  return(download_file_cached(url, dest_file, overwrite = overwrite))
}

#' Download HGNC complete set
#'
#' @param cache_dir Directory to cache the downloaded file. If NULL, proposes .cache/ in R project or current dir.
#' @param overwrite Logical, whether to overwrite existing file without prompting (default FALSE)
#' @return Path to the downloaded file
#' @export
download_hgnc_complete <- function(cache_dir = NULL, overwrite = FALSE) {
  base_cache <- get_cache_dir(cache_dir = cache_dir)
  gtex_cache <- file.path(base_cache, "GTEx")
  if (!dir.exists(gtex_cache)) dir.create(gtex_cache, recursive = TRUE)

  url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
  dest_file <- file.path(gtex_cache, basename(url))
  return(download_file_cached(url, dest_file, overwrite = overwrite))
}

#' Read and clean GTEx GENCODE annotation from cache
#'
#' @param cache_dir Cache directory.
#' @return A data frame of cleaned GENCODE genes.
#' @importFrom dplyr filter arrange group_by distinct ungroup .data
#' @keywords internal
read_gtex_gencode_gtf <- function(cache_dir = NULL) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required. Please install it with BiocManager::install('rtracklayer').")
  }

  base_cache <- get_cache_dir(cache_dir = cache_dir)
  gtex_cache <- file.path(base_cache, "GTEx")
  gtf_file <- file.path(gtex_cache, "gencode.v39.GRCh38.genes.gtf")

  if (!file.exists(gtf_file)) {
    stop("GENCODE GTF file not found in cache (", gtf_file, "). Please run download_gtex_gencode() first.")
  }

  message("Reading and cleaning GTF from ", gtf_file, "...")
  gtf <- rtracklayer::import(gtf_file) |>
    as.data.frame() |>
    dplyr::filter(.data$type == "gene") |>
    dplyr::filter(.data$gene_type %in% c("protein_coding", "lncRNA", "miRNA")) |>
    dplyr::arrange(.data$gene_name, .data$seqnames) |>
    dplyr::group_by(.data$gene_name) |>
    dplyr::distinct(.data$gene_name, .keep_all = TRUE) |>
    dplyr::ungroup()

  return(gtf)
}

#' Export GTEx TPM data to ARACNe format
#'
#' @param tissue Tissue name (e.g., "Esophagus_Gastroesophageal_Junction", "Stomach"). Use the same name as in download_gtex_tpm.
#' @param output_dir Directory to save the exported TSV files.
#' @param cache_dir Directory where GTEx data is cached. If NULL, uses standard .cache detection.
#' @param min_tpm Minimum TPM value used for log2 transformation (default 0.009191).
#' @return A list containing the TPM and log2TPM data frames (invisibly).
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr select rename filter distinct mutate across .data
#' @export
export_gtex_aracne <- function(tissue = "Esophagus_Gastroesophageal_Junction", output_dir = ".", cache_dir = NULL, min_tpm = 0.009191) {
  # 1. Load cleaned GTF
  gtf_nodup <- read_gtex_gencode_gtf(cache_dir = cache_dir)

  # 2. Get TPM file path
  base_cache <- get_cache_dir(cache_dir = cache_dir)
  gtex_cache <- file.path(base_cache, "GTEx")
  tissue_lower <- tolower(tissue)
  file_name <- paste0("gene_tpm_v10_", tissue_lower, ".gct.gz")
  tpm_file <- file.path(gtex_cache, file_name)

  if (!file.exists(tpm_file)) {
    stop("TPM file not found in cache (", tpm_file, "). Please run download_gtex_tpm(tissue='", tissue, "') first.")
  }

  # 3. Read TPM data
  message("Reading TPM data from ", tpm_file, "...")
  # GTEx GCT files have 2 header rows to skip
  dd <- readr::read_tsv(tpm_file, skip = 2, show_col_types = FALSE) |>
    dplyr::select(-"Name") |>
    dplyr::rename(geneSymbol = "Description") |>
    dplyr::filter(.data$geneSymbol %in% gtf_nodup$gene_name) |>
    dplyr::distinct(.data$geneSymbol, .keep_all = TRUE)

  # 4. Prepare output dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 5. Save raw TPM (Symbol format)
  out_prefix <- file.path(output_dir, paste0("GTEx10_", tissue_lower))
  raw_file <- paste0(out_prefix, "_tpm_SYMBOL.tsv")
  readr::write_tsv(dd, raw_file)
  message("Saved TPM SYMBOL file to ", raw_file)

  # 6. Calculate log2TPM and save
  log2TPM_func <- function(x) {
    x[x <= 0] <- min_tpm / 2
    log2(x)
  }

  message("Calculating log2TPM...")
  # Process only numeric columns
  dd_log2 <- dd |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), log2TPM_func))

  log2_file <- paste0(out_prefix, "_log2tpm_SYMBOL.tsv")
  readr::write_tsv(dd_log2, log2_file)
  message("Saved log2TPM SYMBOL file to ", log2_file)

  return(invisible(list(tpm = dd, log2tpm = dd_log2)))
}

#' Create a SLURM script for ARACNe analysis
#'
#' @param tissue Tissue name.
#' @param output_dir Directory where the work is performed.
#' @param aracne_jar Path to the aracne.jar file (relative to the output_dir or absolute).
#' @param partition SLURM partition (default "rpci").
#' @param cpus SLURM cpus-per-task (default 32).
#' @param mem SLURM memory (default "64G").
#' @param bootstrap Number of bootstraps (default 500).
#' @return Path to the created SLURM script (invisibly).
#' @export
create_aracne_slurm <- function(tissue = "Esophagus_Gastroesophageal_Junction", output_dir = ".", aracne_jar = "../ARACNe/aracne.jar",
                                partition = "rpci", cpus = 32, mem = "64G", bootstrap = 500) {
  tissue_lower <- tolower(tissue)
  job_name <- paste0("aracne_", tissue_lower)
  infile <- paste0("GTEx10_", tissue_lower, "_log2tpm_SYMBOL.tsv")
  tf_file <- "TF1600_symbol.tsv"
  out_subfolder <- paste0(tissue_lower, "_aracne_results")

  slurm_content <- c(
    "#!/bin/bash",
    "",
    "#SBATCH --clusters=faculty --partition=rpci --qos=rpci",
    paste0("#SBATCH --nodes=1"),
    paste0("#SBATCH --cpus-per-task=", cpus),
    paste0("#SBATCH --job-name=\"", job_name, "\""),
    paste0("#SBATCH --output=%x_%j.out"),
    paste0("#SBATCH --mem=", mem),
    "",
    "# module list",
    "ulimit -s unlimited",
    "module load java/11.0.16",
    "",
    "echo \"Start Time = `date`\"",
    "tic=`date +%s`",
    "",
    paste0("infile=\"", infile, "\""),
    paste0("outdir=\"", out_subfolder, "\""),
    paste0("tf_file=\"", tf_file, "\""),
    paste0("jar_file=\"", aracne_jar, "\""),
    "",
    "mkdir -p $outdir",
    "",
    "## calculate threshold",
    "java -Xmx5G -jar $jar_file -e $infile -o $outdir --tfs $tf_file --pvalue 1E-8 --seed 1 --calculateThreshold",
    "",
    "## bootstrap",
    paste0("for i in {1..", bootstrap, "}"),
    "do",
    paste0("  java -Xmx5G -jar $jar_file -e $infile -o $outdir --tfs $tf_file --pvalue 1E-8 --threads ", cpus, " --seed $i"),
    "done",
    "",
    "## consolidate",
    "java -Xmx5G -jar $jar_file -o $outdir --consolidate",
    "",
    "echo \"End Time = `date`\"",
    "toc=`date +%s`",
    "elapsedTime=`expr $toc - $tic`",
    "echo \"Elapsed Time = $elapsedTime seconds\"",
    ""
  )

  slurm_file <- file.path(output_dir, paste0("slurm_aracne_", tissue_lower, ".sh"))
  writeLines(slurm_content, slurm_file)
  message("Created SLURM script: ", slurm_file)
  return(invisible(slurm_file))
}
