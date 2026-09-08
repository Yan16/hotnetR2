#' Download NCBI feature table for LDAK / JEME workflows
#'
#' Download the NCBI feature table (GRCh37.p13 by default) into
#' the package cache for downstream LDAK and JEME processing. The function
#' returns the local path to the downloaded file.
#'
#' @param cache_dir Character. Directory to cache the downloaded file. If NULL,
#'   the package will propose a `.cache/LDAK` directory in the project or
#'   current working directory.
#' @param url Character. Optional URL pointing to a feature table file. If
#'   NULL the function will use the NCBI GRCh37.p13 default release URL.
#' @param overwrite Logical. If TRUE, overwrite an existing cached file
#'   without prompting. Default is FALSE.
#' @return Character scalar with the filesystem path to the downloaded file.
#' @examples
#' \dontrun{
#' download_NBCI_feature_table()
#' }
#' @seealso download_JEME
#' @family data-download
#' @export
download_NBCI_feature_table <- function(cache_dir = NULL, url = NULL, overwrite = FALSE) {
  # 1. Determine/Create cache directory (using .cache/LDAK)
  ldak_cache <- get_cache_dir(sub_dir = "LDAK", cache_dir = cache_dir)

  # 2. Prepare URL and local file path
  if (is.null(url)) {
    # RefSeq GRCh37.p13 feature table
    message("No URL provided. Defaulting to NCBI GRCh37.p13 feature table...")
    url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.25-RS_2024_09/GCF_000001405.25_GRCh37.p13_feature_table.txt.gz"
  }
  dest_file <- file.path(ldak_cache, basename(url))

  # 3. Download and return
  return(download_file_cached(url, dest_file, overwrite = overwrite))
}

#' Read NCBI feature table and deduplicate
#'
#' Read NCBI feature table and return standardized, deduplicated tibble
#'
#' Read a downloaded NCBI feature table file into a tibble and remove duplicate
#' entries. The function standardizes key columns so downstream code can rely
#' on consistent names and types.
#'
#' The returned tibble includes at minimum the following standardized columns:
#'
#' - `gene_symbol` (character): gene symbol or name extracted from the source.
#' - `CHR` (character): chromosome with any leading `chr` removed and `X`/`Y`
#'    converted to `"23"`/`"24"`.
#' - `START` (integer): start coordinate (TSS / gene start).
#' - `END` (integer): end coordinate (defaults to `START + 1` when no explicit
#'    end is present in the file).
#'
#' Other annotation columns present in the file are preserved when possible.
#' The function deduplicates rows by `gene_symbol` (final tibble contains one
#' row per unique `gene_symbol`).
#'
#' @param file Character scalar path to the downloaded feature table (gz supported).
#' @param keep_classes (Optional) A character vector of classes to keep. If NULL, all classes are retained.
#' @return A tibble with one row per unique gene symbol and standardized
#'   columns including `gene_symbol`, `CHR`, `START`, and `END`.
#' @examples
#' \dontrun{
#' f <- download_NBCI_feature_table()
#' tbl <- NCBI_feature_table(f)
#' }
#' @importFrom readr read_tsv cols
#' @importFrom dplyr select filter rename mutate case_when arrange relocate distinct .data any_of
#' @export
NCBI_feature_table <- function(file, keep_classes = NULL) {
    # The "class" column in the NCBI feature table indicates the type of gene (e.g., "protein_coding", "lncRNA", "miRNA", etc.). If NULL, all classes are retained.
    if (missing(file) || !is.character(file) || length(file) != 1) {
      stop("`file` must be a single path to the downloaded feature table")
    }
    if (!file.exists(file)) stop("file not found: ", file)
    # the first column of this data (# feature) are all "gene", drop it
    # the second column (class) indicate classes of features (e.g., "protein-coding", "lncRNA", etc). We only keep these three classes.
    df <- readr::read_tsv(file, col_types = readr::cols(.default = "c"), progress = FALSE) |>
      dplyr::select(-1) # drop first column of "gene" features;

    # Optionally filter by classes if provided (NULL means keep all)
    if (!is.null(keep_classes)) {
      df <- df |> dplyr::filter(class %in% keep_classes)
    }

    df <- df |> dplyr::select(-dplyr::any_of(c("assembly", "assembly_unit", "seq_type", "product_accession", "related_accession", "non-redundant_refseq")))

    # Standardize key columns and deduplicate
    df |>
      dplyr::rename(gene_symbol = "symbol", CHR = "chromosome", START = "start", END = "end", gene_name = "name") |>
      dplyr::filter(!.data$CHR %in% c("MT", "M")) |>  # exclude mitochondrial genes
      dplyr::mutate(
        CHR = sub("^chr", "", as.character(.data$CHR), ignore.case = TRUE),
        CHR = dplyr::case_when(
          .data$CHR == "X" ~ "23",
          .data$CHR == "Y" ~ "24",
          TRUE ~ .data$CHR
        ) |> as.integer(),
        START = as.integer(.data$START),
        END = as.integer(.data$END)
      ) |>
      dplyr::filter(!is.na(.data$CHR) & !is.na(.data$START) & !is.na(.data$END) & !is.na(.data$gene_symbol)) |>
      dplyr::arrange(.data$CHR, .data$START, .data$END) |>
      dplyr::select(-dplyr::any_of(c("symbol", "chromosome", "chr", "start", "end"))) |>
      dplyr::relocate("gene_symbol", "CHR", "START", "END", "strand", "gene_name", dplyr::everything()) |>
      # final deduplication by gene_symbol
      # TODO: choose the longest one?
      dplyr::distinct(.data$gene_symbol, .keep_all = TRUE)
  }
# ex:
# nbci_file <- download_NBCI_feature_table()
# ncbi <- NCBI_feature_table(nbci_file)  # this is the same as following full path to the file in cache
# ncbi <- NCBI_feature_table(file.path(get_cache_dir(), "LDAK", "GCF_000001405.25_GRCh37.p13_feature_table.txt.gz"))

#' Download and process JEME data
#'
#' @param type Type of data to download: "lasso" or "elasticnet"
#' @param cache_dir Directory to save downloaded zip files; NULL uses the package cache setting.
#' @param out_dir Directory to save processed files. Defaults to ".cache/HiC" safely found within the project root.
#' @param save_summarized Logical, whether to create and save summarized enhancer/promoter RDS and TSV files (default TRUE)
#' @param overwrite Logical, whether to overwrite existing downloaded zip files (default FALSE).
#' @return A list containing the combined JEME dataframe, unique enhancers, and unique promoters
#' @importFrom dplyr select rename group_by summarise across all_of any_of
#' @importFrom tidyr separate
#' @importFrom readr read_csv write_rds write_tsv
#' @importFrom tools file_path_sans_ext
#' @export
download_JEME <- function(type = c("lasso", "elasticnet"), cache_dir = NULL, out_dir = NULL, save_summarized = TRUE, overwrite = FALSE) {
  type <- match.arg(type)

  # Use cache_dir if provided, otherwise it defaults to download_dir logic in old code
  download_dir <- cache_dir

  # 1. Determine/Create cache directory (default to .cache/JEME)
  if (is.null(out_dir)) {
    out_dir <- get_cache_dir(sub_dir = "JEME", cache_dir = download_dir)
  }

  download_dir <- get_cache_dir(cache_dir = download_dir)

  url <- paste0("https://labs.sbpdiscovery.org/centerandlabs/cancercenter/YipLab/jeme/files/encoderoadmap_", type, ".zip")
  zip_file <- file.path(download_dir, paste0("encoderoadmap_", type, ".zip"))
  unzip_dir <- file.path(download_dir, paste0("encoderoadmap_", type))

  sys_type_dir <- file.path(out_dir, type) # e.g., ".../.cache/JEME/lasso"

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!dir.exists(sys_type_dir)) dir.create(sys_type_dir, recursive = TRUE)

  # 2. Download zip file using cached helper
  download_file_cached(url, zip_file, overwrite = overwrite)

  if (!dir.exists(unzip_dir)) {
    message("Unzipping ", zip_file, " to ", unzip_dir)
    utils::unzip(zip_file, exdir = download_dir)
  }

  files <- list.files(unzip_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No csv files found in ", unzip_dir)

  message("Converting ", length(files), " csv files for ", type, " from ", unzip_dir, " to RDS...")
  purrr::walk(list.files(unzip_dir, pattern = "\\.csv$"), ~ {
    out_rds <- file.path(sys_type_dir, paste0(tools::file_path_sans_ext(.x), ".rds"))
    invisible(read_jeme_csv(file.path(unzip_dir, .x), out_rds = out_rds, simplify_col = FALSE))
  })

  message("Loading and combining all tissues using get_jeme()...")
  res <- get_jeme(method = type, simplified = TRUE)

  ret_list <- list(combined = res)

  if (save_summarized) {
    message("Generating unique enhancer and promoter summaries...")

    nodes <- extract_jeme_nodes(res)
    jeme_enhancer <- nodes$enhancers
    jeme_promoter <- nodes$promoters

    enhancer_rds <- file.path(out_dir, paste0("JEME_enhancer_", type, ".rds"))
    promoter_rds <- file.path(out_dir, paste0("JEME_promoter_", type, ".rds"))

    readr::write_rds(jeme_enhancer, enhancer_rds, compress = "gz")
    readr::write_tsv(jeme_enhancer, paste0(enhancer_rds, ".tsv.gz"))
    readr::write_rds(jeme_promoter, promoter_rds, compress = "gz")
    readr::write_tsv(jeme_promoter, paste0(promoter_rds, ".tsv.gz"))

    ret_list$enhancer <- jeme_enhancer
    ret_list$promoter <- jeme_promoter

    message("Saved summarized files to ", out_dir)
  }

  return(ret_list)
}
