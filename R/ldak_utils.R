## Download function moved to R/download_JEME.R as `download_NBCI_feature_table()`

## promoter_range() moved to R/JEME_access.R


#' Create LDAK promoter boundary .loc file
#'
#' This function processes the NCBI feature table to create a gene boundary file
#' compatible with LDAK. The promoter is the strand-aware TSS encoded as the
#' one-base BED interval `[TSS-1,TSS)`. The final promoter range (including
#' flanks) is defined by LDAK using this one-base interval as an anchor.
#' Note that boundaries are created for all features present in the source file,
#' providing a comprehensive reference not limited to JEME or HiC specific regions.
#'
#' @param feature_file NCBI feature table path; NULL searches the LDAK cache.
#' @param cache_dir Directory where the feature table is cached.
#' @param output_file Path to the output .loc file. If NULL, saves to cache_dir.
#' @return Path to the created .loc file
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr mutate select filter distinct group_by arrange ungroup relocate last_col .data
#' @export
create_ldak_promoter_boundary <- function(feature_file = NULL, cache_dir = NULL, output_file = NULL) {
  base_cache <- get_cache_dir(cache_dir = cache_dir)
  ldak_cache <- file.path(base_cache, "LDAK")
  if (is.null(feature_file)) {
    feature_file <- list.files(ldak_cache, pattern = "feature_table.txt.gz$", full.names = TRUE)[1]
  }
  if (is.na(feature_file) || !file.exists(feature_file)) {
    stop("NCBI feature table not found in ", ldak_cache, ". Please run download_ldak_feature_table() first.")
  }

  message("Reading feature table from ", feature_file, "...")
  raw <- readr::read_tsv(feature_file, show_col_types = FALSE)

  # Processing logic from tmp/create_promoter_boundary.R
  message("Processing gene boundaries for LDAK...")
  ldak <- raw |>
    promoter_range(source = "NCBI") |>
    dplyr::select("symbol", "CHR", "START", "END", "strand") |>
    dplyr::filter(!is.na(.data$CHR)) |>
    dplyr::filter(.data$CHR != "MT") |>
    dplyr::distinct() |>
    dplyr::group_by(.data$symbol, .data$CHR) |>
    dplyr::arrange(.data$symbol, .data$CHR, .data$START) |>
    dplyr::filter(dplyr::row_number() == 1) |>
    dplyr::ungroup()

  # Rename X/Y to 23/24 and ensure numeric for LDAK sorting
  ldak$CHR[ldak$CHR == "X"] <- "23"
  ldak$CHR[ldak$CHR == "Y"] <- "24"
  ldak$CHR <- as.numeric(ldak$CHR)

  ldak <- ldak |>
    dplyr::relocate("strand", .after = dplyr::last_col()) |>
    dplyr::arrange(.data$CHR, .data$START, .data$END)

  if (is.null(output_file)) {
    # Generate default output name based on input file
    base_name <- sub("_feature_table.txt.gz$", "", basename(feature_file))
    output_file <- file.path(get_data_refs_dir(), paste0("LDAK_promoter_", base_name, ".loc"))
  }

  message("Saving LDAK .loc file to ", output_file, "...")
  readr::write_tsv(ldak, output_file, col_names = FALSE)

  return(output_file)
}


#' Export extracted enhancers to LDAK .loc format
#'
#' @param nodes A list containing `enhancers` data frame, as returned by
#'   `extract_jeme_nodes()` or `extract_hic_nodes()`.
#' @param output_file Path to the output .loc file.
#' @param source_type Character. Either "jeme" or "hic" to identify the name column.
#' @param cache_dir Directory to cache the output. If NULL, uses .cache/LDAK/.
#' @return ldak_df data frame containing the enhancer coordinates in LDAK format
#' @importFrom readr write_tsv
#' @importFrom dplyr mutate select arrange distinct case_when across .data
#' @importFrom stringr str_remove
#' @export
export_enhancer_ldak <- function(nodes, output_file = NULL, source_type = c("jeme", "hic"), cache_dir = NULL) {
  source_type <- match.arg(source_type)

  if (is.null(output_file)) {
    output_file <- file.path(get_data_refs_dir(cache_dir), "enhancer_ldak.loc")
  }

  df <- nodes$enhancers
  if (is.null(df) || nrow(df) == 0) {
    stop("No enhancers found in the nodes object.")
  }

  # Identify name column
  name_col <- if (source_type == "jeme") "enhancer" else "Interacting_fragment"
  if (!name_col %in% colnames(df)) {
    stop("Expected column '", name_col, "' not found for source_type '", source_type, "'.")
  }

  # Logic similar to create_*_boundary
  ldak_df <- df |>
    dplyr::select(name = dplyr::all_of(name_col), "CHR", "START", "END") |>
    dplyr::mutate(CHR = stringr::str_remove(.data$CHR, "chr")) |>
    dplyr::mutate(CHR = dplyr::case_when(
      .data$CHR == "X" ~ "23",
      .data$CHR == "Y" ~ "24",
      TRUE ~ .data$CHR
    )) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c("CHR", "START", "END")), as.integer)) |>
    dplyr::arrange(.data$CHR, .data$START, .data$END) |>
    dplyr::distinct()

  message("Saving enhancers to LDAK .loc file: ", output_file)
  readr::write_tsv(ldak_df, output_file, col_names = FALSE)

  return(ldak_df)
}

# Examples:
# nodes <- extract_jeme_nodes(jeme_data)
# export_enhancer_ldak(nodes, source_type = "jeme")

#' Deprecated: Export extracted promoters to LDAK .loc format
#'
#' @description This function is deprecated. Please use \code{\link{jeme_promoter_to_loc}} instead.
#' @param nodes A list containing `promoters` data frame, as returned by
#'   `extract_jeme_nodes()` or `extract_hic_nodes()`.
#' @param output_file Path to the output .loc file. We only write the first 4 columns (name, CHR, START, END) for LDAK compatibility.
#' @param cache_dir Directory to cache the output. If NULL, uses .cache/LDAK/.
#' @return ldak_df data frame containing the promoter coordinates in LDAK format.
#' @export
export_promoter_ldak <- function(nodes, output_file = NULL, cache_dir = NULL) {
  .Deprecated("jeme_promoter_to_loc", package = "hotnetR2")

  if (is.null(output_file)) {
    output_file <- file.path(get_data_refs_dir(cache_dir), "promoter_ldak.loc")
  }

  df <- nodes$promoters
  if (is.null(df) || nrow(df) == 0) {
    stop("No promoters found in the nodes object.")
  }

  # Identify name column: prioritize harmonized symbol, then promoter symbol, then ENSG
  name_col <- dplyr::case_when(
    "promoter_harmonized" %in% colnames(df) ~ "promoter_harmonized",
    "promoter" %in% colnames(df) ~ "promoter",
    "Promoter" %in% colnames(df) ~ "Promoter", # For HiC consistency
    "ENSG" %in% colnames(df) ~ "ENSG",
    TRUE ~ colnames(df)[1]
  )

  # Logic similar to export_enhancer_ldak
  ldak_df <- df |>
    dplyr::select(name = dplyr::all_of(name_col), "CHR", "START", "END") |>
    dplyr::mutate(CHR = stringr::str_remove(.data$CHR, "chr")) |>
    dplyr::mutate(CHR = dplyr::case_when(
      .data$CHR == "X" ~ "23",
      .data$CHR == "Y" ~ "24",
      TRUE ~ .data$CHR
    )) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c("CHR", "START", "END")), as.integer)) |>
    dplyr::arrange(.data$CHR, .data$START, .data$END) |>
    dplyr::distinct()
  # Ensure we only write the first four columns required by LDAK (.loc): name, CHR, START, END
  out_df <- ldak_df |> dplyr::select(dplyr::all_of(c("name", "CHR", "START", "END")))

  message("Saving promoters to LDAK .loc file: ", output_file)
  readr::write_tsv(out_df, output_file, col_names = FALSE)

  return(ldak_df)
}

#' Convert JEME promoter nodes to LDAK .loc format
#'
#' @param promoters A data frame of JEME promoter nodes, as returned by `extract_jeme_nodes()$promoters`.
#' @param output_file Path to the output .loc file. If NULL, defaults to JEME subdirectory in data reference.
#' @param annotation_src Character. Either "NCBI" (use NCBI feature table for ranges) or "jeme" (use coordinates in JEME data).
#' @param ncbi_tbl A data frame from `NCBI_feature_table()`, required if `annotation_src == "NCBI"`.
#' @return A data frame of the formatted LDAK boundary data.
#' @importFrom readr write_tsv
#' @importFrom dplyr mutate select arrange distinct case_when across .data all_of inner_join
#' @export
jeme_promoter_to_loc <- function(promoters, output_file = NULL, annotation_src = c("NCBI", "jeme"), ncbi_tbl = NULL) {
  if (is.null(output_file)) {
    output_file <- file.path(get_data_refs_dir(), "jeme_promoters.loc")
  }
  annotation_src <- match.arg(annotation_src)

  if (annotation_src == "NCBI") {
    # Use NCBI-based promoter ranges
    # Match promoter name to NCBI gene symbol
    message("Using NCBI-based promoter ranges for LDAK .loc file.")
    # ncbi_tbl <- NCBI_feature_table() # expects the same feature table used for create_ldak_promoter_boundary()
    promoter_df <- promoters |>
      dplyr::select(promoter, CHR) |>
      dplyr::inner_join(ncbi_tbl, by = c("promoter" = "gene_symbol", "CHR" = "CHR")) |>
      dplyr::mutate(
        TSS = dplyr::if_else(.data$strand == "+", .data$START, .data$END),
        START = .data$TSS - 1L,
        END = .data$TSS
      ) |>
      dplyr::select(name = "promoter", dplyr::all_of(c("CHR", "START", "END", "strand", "GeneID")))
  } else {
    # Use JEME-provided coordinates directly
    message("Using JEME-provided coordinates for LDAK .loc file.")
    if (!"location" %in% colnames(promoters)) stop("Source df must contain 'location' column for jeme source.")
    promoter_df <- promoters |>
      dplyr::mutate(
        START = as.integer(.data$location) - 1L,
        END = as.integer(.data$location)
      ) |>
      dplyr::select(dplyr::all_of(c("promoter", "CHR", "START", "END", "strand")),
                    name = "promoterFull", dplyr::all_of(c("ENSG", "location"))) |>
      dplyr::distinct()
  }

  promoter_df <- promoter_df |> dplyr::arrange(CHR, START, END)

  promoter_df |>
    dplyr::select(dplyr::all_of(c("name", "CHR", "START", "END", "strand"))) |>
    readr::write_tsv(output_file, col_names = FALSE)

  return(promoter_df)
}

## example usage:
# jeme_nodes <- extract_jeme_nodes(jeme_data)
# ncbi_file <- file.path(get_cache_dir(), "LDAK", "GCF_000001405.25_GRCh37.p13_feature_table.txt.gz")
# ncbi <- NCBI_feature_table(ncbi_file)
# jeme_loc <- jeme_promoter_to_loc(jeme_nodes$promoters, annotation_src = "NCBI", ncbi = ncbi)
# dim(jeme_loc) # [1] 14148    13
##  Use location
# jeme_loc2 <- jeme_promoter_to_loc(jeme_nodes$promoters, annotation_src = "jeme")
# dim(jeme_loc2) # [1] 15673    10


#' Summarize LDAK gene-based analysis results
#'
#' This function reads and summarizes the results from an LDAK-GBAT analysis,
#' including the joined REML results and the estimated genetic contributions.
#' It applies Bonferroni correction to the p-values.
#'
#' @param result_dir Path to the directory containing LDAK results (where remls.all and prs.all.pvalues are located).
#' @param output_file Optional path to save the summary as an RDS file.
#' @return A list containing two data frames:
#'   \item{reml}{The joined REML results with an added `bon.adj` column.}
#'   \item{egc}{The estimated genetic contributions (from prs.all.pvalues) with an added `bon.adj` column.}
#' @importFrom readr read_delim write_rds
#' @importFrom dplyr mutate arrange select starts_with
#' @importFrom stringr str_c
#' @export
#' @examples
#' \dontrun{
#' # Assuming results are in "results/EA_META_1k"
#' summary <- summarize_ldak("results/EA_META_1k")
#'
#' # View top significant genes from REML
#' head(summary$reml[summary$reml$bon.adj < 0.05, ])
#'
#' # View top significant genes from EGC
#' head(summary$egc[summary$egc$bon.adj < 0.05, ])
#' }
summarize_ldak <- function(result_dir, output_file = NULL) {
  if (!dir.exists(result_dir)) {
    stop("Result directory does not exist: ", result_dir)
  }

  reml_path <- file.path(result_dir, "remls.all")
  egc_path <- file.path(result_dir, "prs.all.pvalues")

  if (!file.exists(reml_path)) {
    stop("REML results file not found: ", reml_path)
  }

  message("Reading REML results from ", reml_path, "...")
  # Column types based on LDAK remls.all format
  reml <- readr::read_delim(reml_path, col_types = "cciiinnnnnnnnnn", show_col_types = FALSE)

  if ("LRT_P_Perm" %in% colnames(reml)) {
    reml <- reml |>
      dplyr::mutate(bon.adj = stats::p.adjust(.data$LRT_P_Perm, method = "bonferroni"))
  } else {
    warning("Column 'LRT_P_Perm' not found in REML results. Bonferroni correction skipped for REML.")
  }

  egc <- NULL
  if (file.exists(egc_path)) {
    message("Reading EGC results from ", egc_path, "...")
    egc <- readr::read_delim(egc_path, col_names = c("Gene_name", "P"), show_col_types = FALSE) |>
      dplyr::mutate(bon.adj = stats::p.adjust(.data$P, method = "bonferroni"))
  } else {
    warning("EGC results file not found: ", egc_path)
  }

  res <- list(reml = reml, egc = egc)

  if (!is.null(output_file)) {
    message("Saving summary to ", output_file, "...")
    readr::write_rds(res, output_file)
  }

  return(res)
}

#' Summarize LDAK gene-based analysis results based on Configuration
#'
#' This function replicates the behavior of \code{analyze_ldak_results2.R}. It
#' reads the requested YAML configuration, loads the REML results, merges the
#' minimum p-values, calculates FDR, and performs a range-based join to identify
#' which SNP predictors from the GWAS summary overlap with the flanks of
#' significant genes.
#'
#' @param config_file Path to the LDAK parameter YAML configuration file.
#' @param fdr_cutoff Numeric significance threshold for FDR. Default is 0.05.
#' @param out_top_gene Name/path for the output top genes file. Default is "top_gene.txt".
#' @param out_top_gene2 Name/path for the output SNP-level top genes file. Default is "top_gene2.txt".
#' @return A list containing `results` (the complete merged REML results), `top_hits`
#'   (significant genes below the FDR threshold), and `top_gene2` (the joined SNP-level data).
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr select left_join mutate filter arrange inner_join rename join_by
#' @importFrom tidyr separate
#' @importFrom yaml read_yaml
#' @importFrom stats p.adjust
#' @importFrom rlang .data
#' @export
summarize_ldak_results <- function(config_file = "ldak_parameter.yaml",
                                   fdr_cutoff = 0.05,
                                   out_top_gene = "top_gene.txt",
                                   out_top_gene2 = "top_gene2.txt") {
  if (!file.exists(config_file)) {
    stop("Error: configuration file not found: ", config_file)
  }
  config <- yaml::read_yaml(config_file)

  output_dir <- paste0(config$out_prefix, "_", config$flank_char, "_output")
  reml_file <- file.path(output_dir, "remls.all")
  details_file <- file.path(output_dir, "genes.details")
  used_snps_file <- file.path(output_dir, "genes.predictors.used")
  summary_file <- config$summaries_out
  flank_size <- as.numeric(config$flank_num)

  message(">>> Loading LDAK output files from: ", output_dir)

  if (!file.exists(reml_file)) stop("Error: remls.all not found at ", reml_file)
  results <- readr::read_tsv(reml_file, show_col_types = FALSE)

  if (file.exists(details_file)) {
    message(">>> Merging Min_Pvalue from genes.details...")
    details <- readr::read_tsv(details_file, skip = 2, show_col_types = FALSE)
    if (all(c("Gene_Name", "Min_Pvalue") %in% colnames(details))) {
      details_subset <- details |> dplyr::select("Gene_Name", "Min_Pvalue")
      results <- results |> dplyr::left_join(details_subset, by = "Gene_Name")
    } else {
      warning("Required columns 'Gene_Name' or 'Min_Pvalue' not found in genes.details.")
    }
  } else {
    warning("genes.details not found. Skipping minP merge.")
  }

  if ("LRT_P_Perm" %in% colnames(results)) {
    results <- results |>
      dplyr::mutate(FDR = stats::p.adjust(.data$LRT_P_Perm, method = "fdr")) |>
      dplyr::arrange(.data$FDR, .data$LRT_P_Perm)
  } else {
    stop("Column 'LRT_P_Perm' not found in REML results.")
  }

  top_hits <- results |> dplyr::filter(.data$FDR < fdr_cutoff)

  readr::write_tsv(top_hits, out_top_gene)
  message(">>> Saved significant genes (FDR < ", fdr_cutoff, ") to: ", out_top_gene)

  top_gene2_df <- NULL
  if (file.exists(summary_file) && file.exists(used_snps_file)) {
    message(">>> Subsetting GWAS summaries to match used SNPs...")

    summaries <- readr::read_tsv(summary_file, show_col_types = FALSE)
    used_snps <- readr::read_tsv(used_snps_file, col_names = FALSE, show_col_types = FALSE)
    used_snps_list <- used_snps[[1]]

    summaries_used <- summaries |> dplyr::filter(.data$Predictor %in% used_snps_list)

    # Parse Predictor "CHR:POS" into SNP_Chr and SNP_Pos
    summaries_used <- summaries_used |>
      tidyr::separate(.data$Predictor, into = c("SNP_Chr", "SNP_Pos"), sep = ":", remove = FALSE, convert = TRUE)

    results_flank <- results |>
      dplyr::mutate(
        start_flank = .data$Gene_Start - flank_size,
        end_flank = .data$Gene_End + flank_size
      )

    message(">>> Linking used SNPs to gene boundaries...")
    top_gene2_df <- dplyr::inner_join(
      summaries_used,
      results_flank,
      dplyr::join_by("SNP_Chr" == "Gene_Chr", "SNP_Pos" >= "start_flank", "SNP_Pos" <= "end_flank")
    )

    top_gene2_df <- top_gene2_df |>
      dplyr::rename(Window_Start = "start_flank", Window_End = "end_flank")

    readr::write_tsv(top_gene2_df, out_top_gene2)
    message(">>> Saved SNP-to-Gene mapping to: ", out_top_gene2)

  } else {
    message(">>> Required SNP files missing. Skipping ", out_top_gene2, " generation.")
  }

  message("\n--- Summary ---")
  message("Significant Genes: ", nrow(top_hits))
  if (!is.null(top_gene2_df)) message("Total Predictors Linked: ", nrow(top_gene2_df))

  return(invisible(list(results = results, top_hits = top_hits, top_gene2 = top_gene2_df)))
}
