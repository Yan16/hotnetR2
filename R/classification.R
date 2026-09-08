#' Classify enhancers against the configured promoter windows
#'
#' Reproduces the v2/v3 same-chromosome overlap policy using bedtools. Class 2
#' takes precedence when at least one window fully contains the enhancer.
#' The compatibility conversion deliberately treats stored coordinates as
#' inclusive when constructing the bedtools intervals, matching the frozen
#' classifier scripts. It does not reinterpret historical coordinates.
#'
#' @param config Configuration returned by [read_analysis_config()].
#' @param bedtools Path or command name of the bedtools executable.
#' @return Invisibly, the complete classification tibble. Writes the class
#'   audit, class-specific annotation files, and class-1 LDAK input files.
#' @export
classify_enhancers <- function(config, bedtools = "bedtools") {
  paths <- analysis_paths(config)
  annotations <- paths[["ldak_annotations"]]
  universe <- file.path(annotations, "enhancers.all.details.tsv.gz")
  if (!file.exists(universe)) {
    file.copy(file.path(annotations, "enhancers.details.tsv.gz"), universe)
  }
  enhancer <- readr::read_tsv(universe,
                            show_col_types = FALSE)
  promoter <- readr::read_tsv(file.path(annotations, "promoters.details.tsv.gz"),
                            show_col_types = FALSE)
  settings <- config$enhancer_classification
  temp <- withr::local_tempdir(pattern = "hotnetR2-overlap-")
  make_bed <- function(data, flank) {
    data |>
      dplyr::transmute(chr = paste0("chr", .data$CHR),
                       start = pmax(1, .data$START - flank) - 1,
                       end = .data$END + flank, name = .data$name) |>
      dplyr::arrange(.data$chr, .data$start, .data$end, .data$name)
  }
  enhancer_bed <- file.path(temp, "enhancers.bed")
  promoter_bed <- file.path(temp, "promoters.bed")
  pairs_file <- file.path(temp, "pairs.tsv")
  readr::write_tsv(make_bed(enhancer, settings$enhancer_flank_bp), enhancer_bed,
                   col_names = FALSE)
  readr::write_tsv(make_bed(promoter, settings$promoter_flank_bp), promoter_bed,
                   col_names = FALSE)
  processx::run(bedtools, c("intersect", "-sorted", "-a", enhancer_bed,
                           "-b", promoter_bed, "-wa", "-wb"), stdout = pairs_file)
  pairs <- readr::read_tsv(
    pairs_file,
    col_names = c("chr", "start", "end", "name", "p_chr", "p_start", "p_end", "promoter"),
    col_types = "cddccddc", progress = FALSE
  ) |>
    dplyr::mutate(contained = .data$start >= .data$p_start & .data$end <= .data$p_end) |>
    dplyr::distinct(.data$name, .data$promoter, .data$contained) |>
    dplyr::arrange(.data$name, .data$promoter)
  overlap <- pairs |>
    dplyr::summarise(
      overlapping_promoter_count = dplyr::n(),
      overlapping_promoters = paste(.data$promoter, collapse = ";"),
      fully_containing_promoter_count = sum(.data$contained),
      fully_containing_promoters = paste(.data$promoter[.data$contained], collapse = ";"),
      .by = "name"
    )
  # The initial seven annotation columns are the stable source-detail contract.
  classes <- enhancer |>
    dplyr::select("name", "CHR", "START", "END", "source", "interaction_type", "tissues") |>
    dplyr::left_join(overlap, by = "name", relationship = "one-to-one") |>
    tidyr::replace_na(list(overlapping_promoter_count = 0L, overlapping_promoters = "",
                          fully_containing_promoter_count = 0L, fully_containing_promoters = "")) |>
    dplyr::mutate(
      enhancer_class = dplyr::case_when(
        .data$fully_containing_promoter_count > 0 ~ settings$class2_label,
        .data$overlapping_promoter_count > 0 ~ settings$class3_label,
        TRUE ~ settings$class1_label
      ),
      classification_enhancer_start = pmax(1, .data$START - settings$enhancer_flank_bp),
      classification_enhancer_end = .data$END + settings$enhancer_flank_bp
    ) |>
    dplyr::relocate("enhancer_class", .after = "tissues")
  if (identical(settings$promoter_window, "gene_body")) {
    classes <- classes |>
      dplyr::rename(overlapping_gene_node_count = "overlapping_promoter_count",
                    overlapping_gene_nodes = "overlapping_promoters",
                    fully_containing_gene_node_count = "fully_containing_promoter_count",
                    fully_containing_gene_nodes = "fully_containing_promoters")
  }
  readr::write_tsv(classes, file.path(annotations, "enhancer_classification.tsv.gz"))
  readr::write_tsv(pairs, file.path(annotations, "enhancer_overlap_pairs.tsv.gz"))
  for (i in 1:3) {
    class_name <- settings[[paste0("class", i, "_label")]]
    selected <- dplyr::filter(classes, .data$enhancer_class == class_name)
    prefix <- if (i == 1L) "enhancers" else paste0("enhancers_class", i)
    readr::write_tsv(selected, file.path(annotations, paste0(prefix, ".details.tsv.gz")))
    readr::write_tsv(dplyr::select(selected, "name", "CHR", "START", "END"),
                     file.path(annotations, paste0(prefix, ".loc")), col_names = FALSE)
  }
  invisible(classes)
}
