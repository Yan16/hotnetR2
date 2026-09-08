## code to prepare `key_tissues_jeme` dataset goes here

library(dplyr)
library(readxl)
library(stringr)
library(usethis)

# Path to the original excel file
raw_file <- normalizePath("data-raw/KEY_TISSUES_Nv7_jeme.xlsx", mustWork = FALSE)

if (file.exists(raw_file)) {
  message("Processing ", raw_file)
  key_tissues_jeme <- readxl::read_xlsx(raw_file) |>
    dplyr::rename(desc1 = `desc...5`, desc2 = `desc...11`) |>
    dplyr::mutate(nfile = stringr::str_extract(file, "\\d+"))
  
  usethis::use_data(key_tissues_jeme, overwrite = TRUE)
} else {
  stop("Raw excel file not found at ", raw_file)
}
