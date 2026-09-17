# Internal overlap policy: compare against selected JEME tissues even when
# preparing the all-tissue LDAK annotation universe. PP contacts are untouched.
validate_hic_jeme_overlap <- function(config) {
  settings <- config$regulatory$hic$jeme_overlap
  if (is.null(settings)) return(invisible(NULL))
  if (!is.logical(settings$enabled) || length(settings$enabled) != 1L || is.na(settings$enabled)) {
    stop("regulatory.hic.jeme_overlap.enabled must be true or false", call. = FALSE)
  }
  for (field in c("hic_flank_bp", "jeme_flank_bp")) {
    value <- settings[[field]] %||% 1000
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value != floor(value)) {
      stop("regulatory.hic.jeme_overlap.", field, " must be a nonnegative integer", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Classify HiC enhancers against a selected JEME reference
#'
#' Both input coordinate sets follow the workflow's 1-based inclusive contract.
#' BED conversion makes a shared boundary base an overlap; adjacent bases are not.
#' @param po HiC promoter-other table from the source loader.
#' @param jeme Selected JEME table or tissue-indexed list of tables.
#' @param hic_flank,jeme_flank Nonnegative flank lengths in base pairs.
#' @param bedtools Executable used for sorted interval intersection.
#' @return A list containing classification, overlap pairs and reference intervals.
#'   Audit coordinates are expanded, 1-based and inclusive.
#' @noRd
hic_jeme_overlap_tables <- function(po, jeme, hic_flank = 1000, jeme_flank = 1000,
                                    bedtools = "bedtools") {
  if (is.list(jeme) && !is.data.frame(jeme)) jeme <- dplyr::bind_rows(jeme)
  make_regions <- function(data, id, flank) {
    data |>
      dplyr::transmute(
        chr = toupper(sub("^chr", "", as.character(.data$CHR), ignore.case = TRUE)),
        start = pmax(1, .data$START - flank), end = .data$END + flank,
        name = as.character(.data[[id]])
      ) |>
      dplyr::mutate(chr = dplyr::recode(.data$chr, X = "23", Y = "24", M = "MT")) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$chr, .data$start, .data$end, .data$name)
  }
  hic <- make_regions(po, "Interacting_fragment", hic_flank)
  jeme <- make_regions(jeme, "enhancer", jeme_flank)
  pairs <- tibble::tibble(hic_chr = character(), hic_start = double(), hic_end = double(),
                         hic_enhancer = character(), jeme_chr = character(),
                         jeme_start = double(), jeme_end = double(), jeme_enhancer = character())
  if (nrow(hic) && nrow(jeme)) {
    temp <- withr::local_tempdir(pattern = "hic-jeme-overlap-")
    a <- file.path(temp, "hic.bed")
    b <- file.path(temp, "jeme.bed")
    out <- file.path(temp, "pairs.tsv")
    readr::write_tsv(dplyr::mutate(hic, start = .data$start - 1), a, col_names = FALSE)
    readr::write_tsv(dplyr::mutate(jeme, start = .data$start - 1), b, col_names = FALSE)
    processx::run(bedtools, c("intersect", "-sorted", "-a", a, "-b", b, "-wa", "-wb"), stdout = out)
    pairs <- readr::read_tsv(out, col_names = names(pairs), col_types = "cddccddc", progress = FALSE) |>
      dplyr::mutate(hic_start = .data$hic_start + 1, jeme_start = .data$jeme_start + 1)
  }
  counts <- pairs |>
    dplyr::summarise(overlapping_jeme_count = dplyr::n_distinct(.data$jeme_enhancer), .by = "hic_enhancer")
  classes <- hic |>
    dplyr::rename(hic_enhancer = "name") |>
    dplyr::left_join(counts, by = "hic_enhancer") |>
    tidyr::replace_na(list(overlapping_jeme_count = 0L)) |>
    dplyr::mutate(hic_jeme_class = dplyr::if_else(.data$overlapping_jeme_count == 0L,
                                                "hic_jeme_class1", "hic_jeme_class2"),
                  retained_by_jeme_filter = .data$overlapping_jeme_count == 0L,
                  hic_flank_bp = hic_flank, jeme_flank_bp = jeme_flank)
  list(classification = classes, overlaps = pairs, jeme_reference = jeme)
}

filter_hic_jeme_overlap <- function(hic, config, all_tissues) {
  validate_hic_jeme_overlap(config)
  selected <- config$regulatory$jeme
  jeme <- get_jeme(key_column = selected$key_column, value = selected$value,
                   method = selected$method %||% "lasso", simplified = TRUE,
                   cache_dir = config$cache_dir)
  settings <- config$regulatory$hic$jeme_overlap
  audit <- hic_jeme_overlap_tables(hic$PO, jeme, settings$hic_flank_bp %||% 1000,
                                   settings$jeme_flank_bp %||% 1000)
  audit$classification <- dplyr::mutate(
    audit$classification, jeme_key_column = selected$key_column,
    jeme_selection = paste(unlist(selected$value), collapse = ";")
  )
  scope <- if (all_tissues) "annotations" else "networks"
  for (table in names(audit)) {
    path <- annotation_output_file(config, paste0("hic_jeme_", scope, "_", table, ".tsv.gz"))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(audit[[table]], path)
  }
  retained <- dplyr::filter(audit$classification, .data$retained_by_jeme_filter)
  hic$PO <- dplyr::semi_join(hic$PO, retained, by = c("Interacting_fragment" = "hic_enhancer"))
  hic
}
