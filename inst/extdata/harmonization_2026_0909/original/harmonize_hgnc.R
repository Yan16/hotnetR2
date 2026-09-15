#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript harmonize_hgnc.R HGNC_COMPLETE_SET.tsv OUTPUT_PREFIX")
}

input_file <- args[[1L]]
output_prefix <- args[[2L]]

split_hgnc_values <- function(data, column, priority, label) {
  data |>
    select(
      Hhgnc_id, symbol, Hname, Hlocus_group, Hstatus, Hlocation,
      Halias_symbol, Hprev_symbol, Hgene_group, Hensembl_gene_id
    ) |>
    mutate(id = .data[[column]], priority = priority) |>
    separate_rows(id, sep = "\\|") |>
    mutate(id = trimws(id), mapping_source = label) |>
    filter(!is.na(id), id != "")
}

hgnc <- readr::read_tsv(input_file, show_col_types = FALSE, na = c("", ".")) |>
  rename_with(~ paste0("H", .x)) |>
  rename(
    symbol = Hsymbol,
    Hhgnc_id = Hhgnc_id,
    Halias_symbol = Halias_symbol,
    Hprev_symbol = Hprev_symbol,
    Hensembl_gene_id = Hensembl_gene_id
  ) |>
  mutate(
    priority = 1L,
    id = symbol,
    mapping_source = "current_symbol"
  )

current <- hgnc |>
  select(
    id, priority, mapping_source, Hhgnc_id, symbol, Hname, Hlocus_group,
    Hstatus, Hlocation, Halias_symbol, Hprev_symbol, Hgene_group,
    Hensembl_gene_id
  )

previous <- split_hgnc_values(hgnc, "Hprev_symbol", 2L, "previous_symbol")
aliases <- split_hgnc_values(hgnc, "Halias_symbol", 3L, "alias_symbol")

mapping <- bind_rows(current, previous, aliases) |>
  arrange(id, priority, symbol)

readr::write_tsv(mapping, paste0(output_prefix, "_HGNC_APPENDlong_SORTid.tsv"), na = "")
writexl::write_xlsx(mapping, paste0(output_prefix, "_HGNC_APPENDlong_SORTid.xlsx"))

message("Wrote ", nrow(mapping), " long HGNC mappings")