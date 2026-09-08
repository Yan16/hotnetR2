# Analysis-local promoter harmonization shim. This can be removed once the
# all-tissue mapping workflow is provided directly by hotnetR.

normalise_jeme_ensg <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.[0-9]+$", "", x)
}

collapse_harmonization_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
  paste(x, collapse = ";")
}

build_all_tissue_jeme_ensg_map <- function(jeme_all, gtf) {
  required_jeme <- c("ENSG", "promoter")
  missing_jeme <- setdiff(required_jeme, names(jeme_all))
  if (length(missing_jeme)) {
    stop("All-tissue JEME data lack: ", paste(missing_jeme, collapse = ", "), call. = FALSE)
  }

  if (!"ENSG" %in% names(gtf) && "gene_id" %in% names(gtf)) {
    gtf <- dplyr::rename(gtf, ENSG = "gene_id")
  }
  if (!"gene" %in% names(gtf) && "gene_name" %in% names(gtf)) {
    gtf <- dplyr::rename(gtf, gene = "gene_name")
  }
  missing_gtf <- setdiff(c("ENSG", "gene"), names(gtf))
  if (length(missing_gtf)) {
    stop("GENCODE mapping lacks: ", paste(missing_gtf, collapse = ", "), call. = FALSE)
  }

  requested <- jeme_all |>
    dplyr::transmute(
      ensg_key = normalise_jeme_ensg(.data$ENSG),
      input_promoter = trimws(as.character(.data$promoter))
    ) |>
    dplyr::filter(!is.na(.data$ensg_key), nzchar(.data$ensg_key)) |>
    dplyr::summarise(
      input_promoter_count = dplyr::n_distinct(.data$input_promoter),
      input_promoters = collapse_harmonization_values(.data$input_promoter),
      .by = "ensg_key"
    )

  gtf_map <- gtf |>
    dplyr::transmute(
      ensg_key = normalise_jeme_ensg(.data$ENSG),
      gencode_gene = trimws(as.character(.data$gene))
    ) |>
    dplyr::filter(
      !is.na(.data$ensg_key), nzchar(.data$ensg_key),
      !is.na(.data$gencode_gene), nzchar(.data$gencode_gene)
    ) |>
    dplyr::distinct() |>
    dplyr::summarise(
      gencode_gene_count = dplyr::n_distinct(.data$gencode_gene),
      gencode_genes = collapse_harmonization_values(.data$gencode_gene),
      promoter_harmonized = if (dplyr::n_distinct(.data$gencode_gene) == 1L) {
        dplyr::first(.data$gencode_gene)
      } else {
        NA_character_
      },
      .by = "ensg_key"
    )

  requested |>
    dplyr::left_join(gtf_map, by = "ensg_key") |>
    dplyr::mutate(
      gencode_gene_count = dplyr::coalesce(.data$gencode_gene_count, 0L),
      gencode_genes = dplyr::coalesce(.data$gencode_genes, ""),
      harmonization_status = dplyr::case_when(
        .data$gencode_gene_count == 1L ~ "mapped_by_ensg",
        .data$gencode_gene_count > 1L ~ "ambiguous_ensg_fallback_to_input",
        TRUE ~ "unmapped_ensg_fallback_to_input"
      )
    ) |>
    dplyr::arrange(.data$ensg_key)
}

apply_all_tissue_jeme_ensg_map <- function(jeme, mapping) {
  required <- c("ensg_key", "promoter_harmonized", "harmonization_status")
  missing <- setdiff(required, names(mapping))
  if (length(missing)) {
    stop("JEME ENSG mapping lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(mapping$ensg_key)) {
    stop("JEME ENSG mapping contains duplicate normalized ENSG IDs", call. = FALSE)
  }

  jeme |>
    dplyr::mutate(
      .input_order = dplyr::row_number(),
      ensg_key = normalise_jeme_ensg(.data$ENSG),
      promoter_original = as.character(.data$promoter)
    ) |>
    dplyr::left_join(
      dplyr::select(mapping, dplyr::all_of(required)),
      by = "ensg_key"
    ) |>
    dplyr::mutate(
      promoter = dplyr::coalesce(.data$promoter_harmonized, .data$promoter_original),
      harmonization_status = dplyr::coalesce(
        .data$harmonization_status, "ensg_absent_from_all_tissue_map_fallback_to_input"
      )
    ) |>
    dplyr::arrange(.data$.input_order) |>
    dplyr::select(-dplyr::all_of(c(".input_order", "promoter_harmonized")))
}

harmonize_hic_gene_symbols <- function(hic, alias_file) {
  input <- list(
    PO = hic$PO |> dplyr::rename(Regulator = "Interacting_fragment", Target = "Promoter"),
    PP = hic$PP
  )
  output <- harmonize_hic_gtf(
    input, alias_link_file = alias_file, verbose = FALSE, rm_na = FALSE
  )
  data <- list(
    PO = output$PO |>
      dplyr::mutate(Promoter = .data$Target_harmonized) |>
      dplyr::rename(Interacting_fragment = "Regulator") |>
      dplyr::select(-dplyr::all_of(c("Target", "Target_harmonized"))),
    PP = output$PP |>
      dplyr::mutate(
        Promoter1 = .data$Promoter1_harm,
        Promoter2 = .data$Promoter2_harm
      ) |>
      dplyr::select(-dplyr::all_of(c("Promoter1_harm", "Promoter2_harm")))
  )
  summary <- data.frame(
    hic_po_rows = nrow(hic$PO),
    hic_po_changed = sum(output$PO$Target_harmonized != output$PO$Target, na.rm = TRUE),
    hic_pp_rows = nrow(hic$PP),
    hic_pp_changed = sum(
      output$PP$Promoter1_harm != output$PP$Promoter1 |
        output$PP$Promoter2_harm != output$PP$Promoter2,
      na.rm = TRUE
    )
  )
  list(data = data, summary = summary)
}
