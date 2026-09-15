#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript harmonize_jeme.R JEME_BED HGNC_MAPPING_TSV OUTPUT_XLSX")
}

input_file <- args[[1L]]
hgnc_file <- args[[2L]]
output_file <- args[[3L]]

jeme <- readr::read_tsv(input_file, col_names = FALSE, show_col_types = FALSE) |>
  rename(J1ensg = X4, gene = X5) |>
  mutate(
    J1ensg = sub("\\..*$", "", J1ensg),
    id = gene,
    chromosome = X1
  ) |>
  distinct(gene, .keep_all = TRUE) |>
  select(id, J1ensg, gene, chromosome)

hgnc <- readr::read_tsv(hgnc_file, show_col_types = FALSE) |>
  mutate(across(c(id, symbol), as.character))

resolved <- jeme |>
  left_join(hgnc, by = "id", relationship = "many-to-many") |>
  filter(!is.na(symbol)) |>
  mutate(
    hgnc_chromosome = sub("[pq].*$", "", Hlocation),
    match_priority = case_when(
      id == symbol & priority == 1L & hgnc_chromosome == chromosome ~ 1L,
      id != symbol & priority == 2L & hgnc_chromosome == chromosome ~ 2L,
      id != symbol & priority == 3L & hgnc_chromosome == chromosome ~ 3L,
      TRUE ~ 99L
    )
  ) |>
  group_by(id) |>
  arrange(match_priority, symbol, .by_group = TRUE) |>
  slice_head(n = 1L) |>
  ungroup() |>
  select(id, J1ensg, symbol, priority, Hlocation, Hhgnc_id, Hname,
         Hlocus_group, Hstatus, Hgene_group, Hensembl_gene_id)

writexl::write_xlsx(resolved, output_file)
message("Wrote ", nrow(resolved), " JEME mappings to ", output_file)