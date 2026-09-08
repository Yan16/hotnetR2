#' Download and process HiC data
#'
#' @param cache_dir Directory to save downloaded Excel files; NULL uses the package cache setting.
#' @param out_dir Directory to save processed files. Defaults to ".cache/HiC" safely found within the project root.
#' @param overwrite Logical, whether to overwrite existing downloaded files (default FALSE).
#' @return A list containing the P-O and P-P interaction dataframes
#' @importFrom readxl read_xlsx
#' @importFrom tidyr separate
#' @importFrom dplyr mutate rename
#' @importFrom readr write_rds
#' @importFrom tibble rowid_to_column
#' @export
download_HiC <- function(cache_dir = NULL, out_dir = NULL, overwrite = FALSE) {
  download_dir <- cache_dir

  # 1. Determine/Create cache directory (default to .cache/HiC)
  if (is.null(out_dir)) {
    out_dir <- get_cache_dir(sub_dir = "HiC", cache_dir = download_dir)
  }

  download_dir <- get_cache_dir(cache_dir = download_dir)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # URLs from Jung et al. (Nature Genetics 2019) Supplementary Tables
  url_s3 <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0494-8/MediaObjects/41588_2019_494_MOESM3_ESM.xlsx"
  url_s4 <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0494-8/MediaObjects/41588_2019_494_MOESM4_ESM.xlsx"

  file_s3 <- file.path(download_dir, "41588_2019_494_MOESM3_ESM.xlsx")
  file_s4 <- file.path(download_dir, "41588_2019_494_MOESM4_ESM.xlsx")

  # 2. Download files using cached helper
  download_file_cached(url_s3, file_s3, overwrite = overwrite)
  download_file_cached(url_s4, file_s4, overwrite = overwrite)

  message("Processing Table S3 (P-O interactions)...")
  s3 <- readxl::read_xlsx(file_s3, skip = 1) |>
    tidyr::separate("Interacting_fragment", into = c("CHR", "START", "END"), remove = FALSE) |>
    dplyr::mutate(START = as.integer(.data$START), END = as.integer(.data$END))

  out_s3 <- file.path(out_dir, "Jung_HiC_P-O.rds")
  readr::write_rds(s3, out_s3, compress = "gz")
  message("Saved P-O interactions to ", out_s3)

  message("Processing Table S4 (P-P interactions)...")
  s4 <- readxl::read_xlsx(file_s4, skip = 1) |>
    tibble::rowid_to_column("row_id") |>
    dplyr::rename(Promoter1 = "Promoter...1", Promoter2 = "Promoter...2")

  out_s4 <- file.path(out_dir, "Jung_HiC_P-P.rds")
  readr::write_rds(s4, out_s4, compress = "gz")
  message("Saved P-P interactions to ", out_s4)

  return(list(`P-O` = s3, `P-P` = s4))
}
