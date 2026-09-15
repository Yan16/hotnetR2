#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: Rscript harmonize_hic.R HIC_PO_XLSX HIC_PP_XLSX HGNC_MAPPING_TSV OUTPUT_PREFIX")
}

po_file <- args[[1L]]
pp_file <- args[[2L]]
hgnc_file <- args[[3L]]
output_prefix <- args[[4L]]

typo_map <- c(
  `1-Dec` = "DELEC1", `1-Mar` = "MARCHF1", `1-Sep` = "SEPTIN1",
  `2-Mar` = "MARCHF2", `3-Mar` = "MARCHF3", `3-Sep` = "SEPTIN3",
  `4-Mar` = "MARCHF4", `6-Mar` = "MARCHF6", `6-Sep` = "SEPTIN6",
  `7-Mar` = "MARCHF7", `7-Sep` = "SEPTIN7", `8-Mar` = "MARCHF8",
  `9-Sep` = "SEPTIN9", `10-Mar` = "MARCHF10", `11-Mar` = "MARCHF11",
  `11-Sep` = "SEPTIN11"
)

repair_symbols <- function(x) {
  x <- gsub(" ", "", x, fixed = TRUE)
  repaired <- unname(typo_map[x])
  ifelse(!is.na(repaired), repaired, x)
}

hgnc <- readr::read_tsv(hgnc_file, show_col_types = FALSE) |>
  mutate(id = as.character(id), symbol = as.character(symbol))

po <- readxl::read_excel(po_file) |>
  mutate(Pro = repair_symbols(Pro)) |>
  distinct(Pro, .keep_all = TRUE) |>
  select(Pro) |>
  rename(id = Pro) |>
  left_join(hgnc |> select(id, symbol, priority, Hlocation, Hhgnc_id, Hname,
                          Hlocus_group, Hstatus, Hgene_group, Hensembl_gene_id),
            by = "id", relationship = "many-to-many") |>
  group_by(id) |>
  arrange(priority, symbol, .by_group = TRUE) |>
  slice_head(n = 1L) |>
  ungroup() |>
  select(id, symbol, priority, Hlocation, Hhgnc_id, Hname, Hlocus_group,
         Hstatus, Hgene_group, Hensembl_gene_id)

pp_input <- readxl::read_excel(pp_file)
if (!"Pro" %in% names(pp_input)) {
  if (!"Promoter" %in% names(pp_input)) {
    stop("HiC-PP input must contain either `Pro` or `Promoter`")
  }
  pp_input <- pp_input |> rename(Pro = Promoter)
}

pp <- pp_input |>
  select(Tissue_type, Pro, B) |>
  tidyr::separate_rows(Pro, sep = ";") |>
  tidyr::separate_rows(B, sep = ";") |>
  mutate(
    Pro = repair_symbols(Pro),
    B = repair_symbols(B)
  ) |>
  left_join(hgnc |> select(id, symbol), by = c("Pro" = "id"), relationship = "many-to-many") |>
  rename(Pro_sym = symbol) |>
  left_join(hgnc |> select(id, symbol), by = c("B" = "id"), relationship = "many-to-many") |>
  rename(B_sym = symbol) |>
  group_by(Tissue_type, Pro, B) |>
  summarise(Pro_sym = first(na.omit(Pro_sym), default = Pro),
            B_sym = first(na.omit(B_sym), default = B), .groups = "drop")

writexl::write_xlsx(po, paste0(output_prefix, "_hiC_PO.xlsx"))
writexl::write_xlsx(pp, paste0(output_prefix, "_hiC_PP_pairs.xlsx"))
message("Wrote HiC-PO and HiC-PP outputs with prefix ", output_prefix)