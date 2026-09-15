# External data preparation, not loaded into the hotnetR2 namespace.
# Source this file and call build_alias_data(); see the adjacent README.
build_alias_data <- function(source_dir, legacy_file, output_dir) {
  jeme <- readxl::read_excel(file.path(source_dir, "final_harmonize_jemeF.xlsx"))
  po <- readxl::read_excel(file.path(source_dir, "final_harmonize_hiC_PO.xlsx"))
  pp <- readxl::read_excel(file.path(source_dir, "final_harmonize_hiC_PP_pairs.xlsx"))
  observations <- dplyr::bind_rows(
    dplyr::transmute(jeme, gene = .data$id, Hsym = .data$symbol, source = "JEME"),
    dplyr::transmute(po, gene = .data$id, Hsym = .data$symbol, source = "HiC_PO"),
    dplyr::transmute(pp, gene = .data$Pro, Hsym = .data$Pro_sym, source = "HiC_PP_Pro"),
    dplyr::transmute(pp, gene = .data$B, Hsym = .data$B_sym, source = "HiC_PP_B")
  ) |>
    dplyr::mutate(dplyr::across(c("gene", "Hsym"), stringr::str_trim))
  if (anyNA(observations$gene) || anyNA(observations$Hsym) ||
      any(!nzchar(observations$gene)) || any(!nzchar(observations$Hsym))) {
    stop("Curated workbook mappings contain missing labels")
  }
  overrides <- tibble::tribble(~gene, ~selected_symbol,
    "C11orf48", "LBHD1", "MEGT1", "LY6G6D")
  resolved <- observations |>
    dplyr::left_join(overrides, by = "gene", relationship = "many-to-one") |>
    dplyr::mutate(original_symbol = .data$Hsym,
      Hsym = dplyr::coalesce(.data$selected_symbol, .data$Hsym))
  conflicts <- resolved |>
    dplyr::distinct(.data$gene, .data$Hsym) |>
    dplyr::count(.data$gene) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(conflicts)) stop("Unresolved curated aliases: ", paste(conflicts$gene, collapse = ", "))
  # Metadata comes from the selected final symbol, never the rejected alternative.
  metadata <- dplyr::bind_rows(jeme, po) |>
    dplyr::select("symbol", "Hhgnc_id", "Hname", "Hlocation", "Hlocus_group") |>
    dplyr::distinct() |>
    dplyr::arrange(.data$symbol, .data$Hhgnc_id) |>
    dplyr::distinct(.data$symbol, .keep_all = TRUE)
  mapping <- resolved |>
    dplyr::summarise(source = paste(sort(unique(.data$source)), collapse = ";"),
                     .by = c("gene", "Hsym")) |>
    dplyr::left_join(metadata, by = c("Hsym" = "symbol"), relationship = "many-to-one") |>
    dplyr::arrange(.data$gene)
  legacy <- tibble::as_tibble(readxl::read_excel(legacy_file))
  # Keep baseline rows and ENSG information for identifiers outside the new keys.
  # Replace target HGNC metadata on covered rows; unavailable Entrez IDs are NA.
  covered <- legacy |>
    dplyr::semi_join(mapping, by = "gene") |>
    dplyr::select(-dplyr::all_of(c("Hsym", "Hhgnc_id", "Hname", "Hlocation", "Hty", "Hentrez_id"))) |>
    dplyr::left_join(mapping, by = "gene", relationship = "many-to-one") |>
    dplyr::mutate(Hty = .data$Hlocus_group, Hentrez_id = NA_real_) |>
    dplyr::select(dplyr::all_of(names(legacy)))
  added <- mapping |>
    dplyr::anti_join(legacy, by = "gene") |>
    dplyr::transmute(gene = .data$gene, Hsym = .data$Hsym,
      Hhgnc_id = .data$Hhgnc_id, Hname = .data$Hname, Hlocation = .data$Hlocation,
      Hty = .data$Hlocus_group)
  alias <- dplyr::bind_rows(dplyr::anti_join(legacy, mapping, by = "gene"), covered, added) |>
    dplyr::arrange(.data$gene, .data$ensg)
  alias_nodup <- dplyr::distinct(alias, .data$gene, .keep_all = TRUE)
  # The existing data file names intentionally load objects alias / alias_nodup.
  for (object in c("alias", "alias_nodup")) {
    value <- get(object)
    attr(value, "mapping_direction") <- "raw_to_approved"
    attr(value, "harmonization_release") <- "2026-09-09+curation-2026-09-14"
    assign(object, value)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  save(alias, file = file.path(output_dir, "alias_link.rda"), compress = "xz", version = 2)
  save(alias_nodup, file = file.path(output_dir, "alias_link_nodup.rda"), compress = "xz", version = 2)
  saveRDS(alias, file.path(output_dir, "alias_link.rds"), version = 2)
  saveRDS(alias_nodup, file.path(output_dir, "alias_link_nodup.rds"), version = 2)
  readr::write_tsv(mapping, file.path(output_dir, "curated_mapping.tsv.gz"))
  readr::write_tsv(dplyr::distinct(dplyr::filter(resolved, !is.na(.data$selected_symbol)),
    .data$gene, .data$source, .data$original_symbol, .data$Hsym),
    file.path(output_dir, "manual_resolution_audit.tsv"))
  readr::write_tsv(tibble::tibble(metric = c("legacy_rows", "curated_keys", "new_keys", "alias_rows", "alias_nodup_rows"),
    count = c(nrow(legacy), nrow(mapping), nrow(added), nrow(alias), nrow(alias_nodup))),
    file.path(output_dir, "build_summary.tsv"))
  inputs <- c(list.files(source_dir, pattern = "[.]xlsx$", full.names = TRUE), legacy_file)
  readr::write_tsv(tibble::tibble(file = basename(inputs),
    sha256 = purrr::map_chr(inputs, ~ digest::digest(file = .x, algo = "sha256"))),
    file.path(output_dir, "source_manifest.tsv"))
  invisible(list(alias = alias, alias_nodup = alias_nodup, mapping = mapping))
}
