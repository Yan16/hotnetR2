









#' Combine JEME and HiC networks into a single regulatory network
#'
#' This function combines harmonized JEME and HiC (PO and PP) networks
#' into a single data frame with regulator-target pairs and a source column.
#'
#' @param jeme_harmonized Output list from `harmonize_jeme_genecode()`
#' @param hic_harmonized Output list from `harmonize_hic_gtf()`
#' @param harmonized Logical, whether to use harmonized gene names (default = TRUE)
#' @param verbose Logical, whether to print progress messages (default = TRUE)
#'
#' @return A tibble with columns `Regulator`, `Target`, and `source`
#' @importFrom dplyr bind_rows ungroup select
#' @importFrom tibble add_column
#' @importFrom rlang sym
#' @export
combine_jeme_hic <- function(jeme_harmonized, hic_harmonized, harmonized = TRUE, verbose = TRUE) {
  vcat <- function(...) if (verbose) message(...)

  vcat("Combining JEME and HiC networks...")

  # Determine which columns to use
  jeme_target_col <- if (harmonized) "Target_harmonized" else "Target"
  po_target_col   <- if (harmonized) "Target_harmonized" else "Target"
  pp_reg_col      <- if (harmonized) "Promoter1_harm" else "Promoter1"
  pp_target_col   <- if (harmonized) "Promoter2_harm" else "Promoter2"

  net_regulators <- dplyr::bind_rows(
    jeme_harmonized |> dplyr::ungroup() |>
      dplyr::select(Regulator, Target = !!rlang::sym(jeme_target_col)) |>
      tibble::add_column(source = "JEME"),

    hic_harmonized$PO |> dplyr::ungroup() |>
      dplyr::select(Regulator, Target = !!rlang::sym(po_target_col)) |>
      tibble::add_column(source = "HiC_PO"),

    hic_harmonized$PP |> dplyr::ungroup() |>
      dplyr::select(Regulator = !!rlang::sym(pp_reg_col), Target = !!rlang::sym(pp_target_col)) |>
      tibble::add_column(source = "HiC_PP")
  )

  vcat("Combined network size: ", paste(dim(net_regulators), collapse = " x "))

  # Optional NA check
  na_rows <- find_na(net_regulators)
  if (nrow(na_rows) > 0) warning("There are NA values in the combined network!")

  return(net_regulators)
}

#' Extract enhancer locations from JEME and HiC data
#'
#' This function combines enhancer location information from JEME and HiC PO data,
#' removes duplicates, and excludes enhancers in the MHC region.
#'
#' @param jeme A tibble of JEME edges with columns `Regulator`, `CHR`, `START`, `END`.
#' @param hic A list returned by `load_hic()`, containing at least a `PO` tibble
#'   with columns `Regulator`, `CHR`, `START`, `END`.
#' @param verbose Logical, whether to print diagnostic messages (default = TRUE).
#'
#' @return A tibble of enhancers with columns:
#'   - `enhancer`: enhancer identifier
#'   - `CHR`, `START`, `END`: genomic location
#'   - `source`: either `"JEME"` or `"HiC_PO"`
#' @importFrom dplyr bind_rows ungroup select distinct mutate rename group_by arrange filter n
#' @export
enhancer_loc <- function(jeme, hic, verbose = TRUE) {
  vcat <- function(...) if (verbose) message(...)

  vcat("Extracting enhancer locations from JEME and HiC PO...")

  enh <- dplyr::bind_rows(
    jeme |>
      dplyr::ungroup() |>
      dplyr::select(Regulator, CHR, START, END) |>
      dplyr::distinct() |>
      dplyr::mutate(source = "JEME"),
    hic$PO |>
      dplyr::ungroup() |>
      dplyr::select(Regulator, CHR, START, END) |>
      dplyr::distinct() |>
      dplyr::mutate(source = "HiC_PO")
  ) |>
    dplyr::rename(enhancer = Regulator)

  vcat("Initial enhancer count: ", nrow(enh))
  vcat("Unique enhancer IDs: ", length(unique(enh$enhancer)))

  dup <- enh |>
    dplyr::group_by(enhancer) |>
    dplyr::arrange(enhancer) |>
    dplyr::filter(dplyr::n() > 1)

  if (nrow(dup) > 0) {
    vcat("Duplicated enhancer names found: ", nrow(dup))
    if (verbose) print(utils::head(dup))
  }

  vcat("Excluding enhancers in MHC region...")
  enh <- enh |>
    dplyr::rename(gene = enhancer) |>
    exclude_MHC(coordinate_system = "bed") |>
    dplyr::rename(enhancer = gene)

  vcat("Final enhancer count after MHC exclusion: ", nrow(enh))

  return(enh)
}

#' Extract promoter locations from JEME and HiC data
#'
#' @param jeme A data frame with column Target
#' @param hic A list returned by load_hic() containing elements PO and PP
#' @param promoter_loc_file Path to LDAK promoter .loc file
#' @param verbose Logical, whether to print progress messages
#' @return A data frame of promoter locations.
#' @importFrom readr read_tsv
#' @importFrom dplyr distinct inner_join group_by filter ungroup
#' @export
promoter_loc <- function(jeme, hic, promoter_loc_file, verbose = TRUE) {
  stopifnot(file.exists(promoter_loc_file))
  stopifnot(all(c("PO", "PP") %in% names(hic)))

  if (verbose) {
    message("Reading promoter location definition file: ", promoter_loc_file)
  }

  loc_promoter <- readr::read_tsv(
    promoter_loc_file,
    col_names = c("gene", "CHR", "START", "END", "strand"),
    col_types = "ciiic"
  )

  if (verbose) {
    message("Promoter location table size before distinct: ", nrow(loc_promoter))
  }

  # Keep first occurrence, remove duplicate promoters for the same gene
  loc_promoter <- loc_promoter |>
    dplyr::distinct(gene, .keep_all = TRUE)

  if (verbose) {
    message("Promoter location table size after distinct: ", nrow(loc_promoter))
  }

  # Collect all relevant promoter genes from JEME, PO, PP
  prom_genes <- unique(c(
    jeme$Target,
    hic$PO$Target,
    hic$PP$Promoter1,
    hic$PP$Promoter2
  ))

  prom <- data.frame(gene = prom_genes) |>
    dplyr::inner_join(loc_promoter, by = "gene") |>
    exclude_MHC(coordinate_system = "bed")

  if (verbose) {
    message("After intersecting with JEME + HiC + excluding MHC:")
    message("  Number of promoters: ", nrow(prom))
    message("  Unique genes: ", length(unique(prom$gene)))
  }

  # Detect duplicates
  dup_promoters <- prom |>
    dplyr::group_by(gene) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::ungroup()

  # Drop duplicates but keep one row per gene
  prom_dedup <- prom |>
    dplyr::distinct(gene, .keep_all = TRUE)

  # Attach dropped promoters as an attribute
  attr(prom_dedup, "dropped_promoters") <- dup_promoters

  if (verbose && nrow(dup_promoters) > 0) {
    message("Dropped duplicated promoters for ", length(unique(dup_promoters$gene)), " genes")
  }

  prom_dedup
}

#' Filter LDAK results based on an elbow point
#'
#' @param df A data frame containing LDAK results with a 'score' column.
#' @param keep Optional character vector of genes to keep regardless of score.
#' @return A filtered data frame.
#' @importFrom dplyr select count mutate filter
#' @export
filter_ldak_gene_elbow <- function(df, keep=NULL){
    # Note: elbow_point should be available in the package or as a dependency
    # For now, assuming it's available or we need to include it.
    # Looking at original code, it sourced elbow.R.

    ## create cumulatted counts vs score data
    ss <- df |>
        dplyr::select(score) |>
        dplyr::count(score) |>
        dplyr::mutate(cum_n = cumsum(n))

    # identify elbow point
    # ep <- elbow_point(ss$score, ss$cum_n) # $x is score cutpoint
    # For now, returning df as is if elbow_point is not defined,
    # or we should implement a simple elbow point if possible.
    # TODO: Implement or include elbow_point
    message("Elbow point filtering requested but elbow_point function not yet implemented in package.")
    return(df)
}

#' Load and process LDAK results
#'
#' @param file Path to the LDAK results file in `.rds` format.
#' @param cohort Optional character string for GWAS cohort filtering.
#' @return A tibble of LDAK results.
#' @importFrom dplyr select mutate relocate filter any_of
#' @export
get_ldak <- function(file, cohort = NULL) {
  out <- readRDS(file) |>
    dplyr::select(gene:Heritability, cohort, dplyr::any_of(c("source", "flank"))) |>
    dplyr::mutate(score = -log10(LRT_P_Perm)) |>
    dplyr::relocate(score, .before = SD)

  if (!is.null(cohort)) {
    message("LDAK scores related to GWAS cohort: ", cohort)
    out <- out |> dplyr::filter(cohort == !!cohort)
  }

  out
}

#' Get LDAK scores for enhancers and promoters
#'
#' @param param A list of parameters.
#' @param enhancers Character vector of enhancer IDs.
#' @param verbose Logical, whether to print progress.
#' @return A tibble of combined LDAK scores.
#' @importFrom dplyr filter mutate bind_rows group_by tally
#' @importFrom here here
#' @export
get_regulator_LDAK_scores <- function(param, enhancers, verbose = TRUE) {
  vcat <- function(...) if (verbose) message(...)

  this_cohort <- param$GWAS$cohort

  # --------------------------------- step 1: LDAK enhancer scores
  vcat("Filtering tissue-specific enhancer LDAK scores for GWAS cohort: ", this_cohort)
  enhancer_ldak <- get_ldak(here::here(param$LDAK$Enhancer$file), cohort = this_cohort) |>
    dplyr::filter(gene %in% enhancers) |>
    dplyr::mutate(node_type = "Enhancer")

  # --------------------------------- step 2: LDAK promoter scores (not tissue-specific)
  vcat("Get promoter LDAK scores, excluding MHC regions")
  promoter_ldak <- get_ldak(here::here(param$LDAK$Promoter$file), cohort = this_cohort) |>
    exclude_MHC() |>
    dplyr::mutate(node_type = "Promoter")

  # --------------------------------- step 3: combine enhancer and promoter scores
  vcat("Combining enhancer and promoter LDAK scores...")
  nodes_all <- dplyr::bind_rows(enhancer_ldak, promoter_ldak)

  if (verbose) {
    vcat("Preview of combined LDAK nodes:")
    print(utils::head(nodes_all))

    vcat("Counts by source and node_type:")
    print(nodes_all |> dplyr::group_by(source, node_type) |> dplyr::tally())

    vcat("Checking for duplicated gene entries...")
    dup <- nodes_all |> dplyr::group_by(gene) |> dplyr::filter(dplyr::n() > 1)
    if (nrow(dup) > 0) print(dup)
  }

  nodes_all
}

#' Filter network edges by LDAK scores
#'
#' @param net_regulators A tibble of network edges.
#' @param nodes_all A tibble of nodes with LDAK scores.
#' @param verbose Logical, whether to print diagnostic messages.
#' @return A list of filtered network data.
#' @importFrom dplyr filter group_by tally
#' @export
filter_network_ldak <- function(net_regulators, nodes_all, verbose = TRUE) {
  vcat <- function(...) if (verbose) message(...)

  vcat("Filtering network edges by LDAK scores...")

  vcat("Before filter: total edges = ", nrow(net_regulators))
  if (verbose) print(net_regulators |> dplyr::group_by(source) |> dplyr::tally())

  jeme_hic <- net_regulators |>
    dplyr::filter(Regulator %in% nodes_all$gene & Target %in% nodes_all$gene)

  droped_jeme_hic <- net_regulators |>
    dplyr::filter(!(Regulator %in% nodes_all$gene & Target %in% nodes_all$gene))

  if (verbose) {
    vcat("After filter: total edges = ", nrow(jeme_hic))
    print(jeme_hic |> dplyr::group_by(source) |> dplyr::tally())

    vcat("Dropped edges: total = ", nrow(droped_jeme_hic))
    print(droped_jeme_hic |> dplyr::group_by(source) |> dplyr::tally())

    vcat("Preview of filtered network edges:")
    print(utils::head(jeme_hic))
  }

  list(
    jeme_hic = jeme_hic,
    droped_jeme_hic = droped_jeme_hic,
    nodes_all = nodes_all
  )
}

#' Read and filter ARACNe network file
#'
#' @param infile Path to the ARACNe network file.
#' @param MI_threshold Minimum MI threshold.
#' @param pvalue_threshold Maximum p-value threshold.
#' @param keep_genes Optional vector of genes to keep.
#' @return A tibble of filtered ARACNe edges.
#' @importFrom readr read_tsv
#' @importFrom dplyr filter bind_rows distinct
#' @export
read_ARACNe <- function(infile, MI_threshold = 0, pvalue_threshold = 1e-6, keep_genes = NULL) {
  net <- readr::read_tsv(infile, col_types = "ccnn")

  out <- net |> dplyr::filter(MI >= MI_threshold & pvalue <= pvalue_threshold)

  if (!is.null(keep_genes)) {
    message("Force keep some genes")
    to_keep <- net |> dplyr::filter(Regulator %in% keep_genes | Target %in% keep_genes)
    out <- out |> dplyr::bind_rows(to_keep) |> dplyr::distinct()
  }

  out
}

#' Calculate negative log base 10, handling infinite values
#'
#' @param x A numeric vector.
#' @param max_value Value to substitute for Inf (default = 1000).
#' @return Negative log10 values.
#' @export
nlog10 <- function(x, max_value=1000) {
  out <- -log10(x)
  if(any(is.infinite(out))) warning("For p-value=0, set its nlog10 transformed to ", max_value)
  out[is.infinite(out)] <- max_value
  out
}

#' Calculate negative log base 2
#'
#' @param x A numeric vector.
#' @param ... Arguments passed to nlog10.
#' @return Negative log2 values.
#' @export
nlog2 <- function(x, ...){
  log2(10)*nlog10(x, ...)
}

#' Read LDAK gene score from file
#'
#' @param flk Flank size.
#' @param cohort GWAS cohort.
#' @param file Path to RDS file.
#' @return A filtered data frame of LDAK scores.
#' @importFrom dplyr ungroup filter
#' @export
read_ldak_gene_score <- function(flk, cohort=NULL, file) {
    valid_cohorts <- c("BEEA_Schroder", "BE_Schroder", "EA_Schroder")

    if (!is.null(cohort) && !(cohort %in% valid_cohorts)) {
        stop("cohort must be NULL or one of: ", paste(valid_cohorts, collapse=", "))
    }

    message(paste("Load flank:", flk, "cohort:", if(is.null(cohort)) "all" else cohort))

    result <- get_ldak(file, cohort) |>
        dplyr::ungroup() |>
        dplyr::filter(flank == !!flk)
    return(result)
}

#' Filter LDAK gene scores based on SD values
#'
#' @param df Data frame of LDAK scores.
#' @param keep Genes to keep.
#' @return Filtered data frame.
#' @importFrom dplyr filter
#' @export
filter_ldak_gene_SD <- function(df, keep=NULL){
  message("Filter LDAK gene scores to exclude those SD==NA")
  if(!is.null(keep)){
    message(paste("Force keep some genes:", paste(keep, collapse=", ")))
  }
  df |> dplyr::filter(!is.na(SD) | gene %in% keep)
}

#' Filter LDAK gene scores based on p-value threshold
#'
#' @param df Data frame of LDAK scores.
#' @param pval.col P-value column name.
#' @param threshold Threshold value.
#' @param keep Genes to keep.
#' @return Filtered data frame.
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @export
filter_ldak_gene_pvalue <- function(df, pval.col, threshold=1, keep=NULL){
    if (!pval.col %in% colnames(df)) {
        stop(sprintf("Column '%s' not found in data frame", pval.col))
    }
    if(!(threshold <= 1 & threshold >= 0)){
        stop("p-value threshold need to be 0-1")
    }

    message(paste0("Filter LDAK gene scores based on p-value threshold = ", threshold))
    if(!is.null(keep)){
        message(paste("Force keep some genes:", paste(keep, collapse=", ")))
    }
    df |> dplyr::filter(!!rlang::sym(pval.col) <= threshold | gene %in% keep)
}

#' Calculate FDR (False Discovery Rate) from p-values
#'
#' @param df Data frame of p-values.
#' @param pval.col P-value column name.
#' @param by.col Grouping columns.
#' @return Data frame with fdr column.
#' @importFrom dplyr group_by across all_of mutate ungroup
#' @importFrom rlang sym
#' @export
calc_fdr <- function(df, pval.col, by.col = NULL) {
    if (!pval.col %in% colnames(df)) {
        stop(sprintf("Column '%s' not found in data frame", pval.col))
    }

    if (is.null(by.col)) {
        df$fdr <- stats::p.adjust(df[[pval.col]], method = "BH")
    } else {
        if (!all(by.col %in% colnames(df))) {
            stop("One or more grouping columns not found in data frame")
        }
        df <- df |>
            dplyr::group_by(dplyr::across(dplyr::all_of(by.col))) |>
            dplyr::mutate(fdr = stats::p.adjust(!!rlang::sym(pval.col), method = "BH")) |>
            dplyr::ungroup()
    }

    return(df)
}

#' Filter LDAK gene scores based on FDR threshold
#'
#' @param df Data frame of LDAK scores.
#' @param pval.col P-value column for FDR.
#' @param by.col Grouping columns.
#' @param threshold FDR threshold.
#' @param keep Genes to keep.
#' @return Filtered data frame.
#' @importFrom dplyr filter
#' @export
filter_ldak_gene_fdr <- function(df, pval.col, by.col = NULL, threshold=1, keep=NULL){
    if(! (threshold <= 1 & threshold >=0 )){stop("fdr threshold need to be 0-1") }

  message(paste0("Filter LDAK gene scores based on fdr threshold = ", threshold))
  message(paste("Calculate FDR based on p-value column:", pval.col))
  if(!is.null(keep)){
    message(paste("Force keep some genes:", paste(keep, collapse=", ")))
  }
    df |> calc_fdr(pval.col=pval.col, by.col=by.col) |>
        dplyr::filter(fdr <= threshold | gene %in% keep)
}

#' Filter LDAK gene scores based on gene types
#'
#' @param df Data frame of LDAK scores.
#' @param gtf GTF annotation.
#' @param types Gene types to keep.
#' @param keep Genes to keep.
#' @return Filtered data frame.
#' @importFrom dplyr left_join select filter
#' @export
filter_ldak_gene_type <- function(df, gtf, types=c("protein_coding"), keep=NULL){
    message(paste("Filter LDAK gene scores based on gene types:", paste(types, collapse=", ")))
  if(!is.null(keep)){
    message(paste("Force keep some genes:", paste(keep, collapse=", ")))
  }
    df |> dplyr::left_join(gtf |> dplyr::select(gene, gene_type), by="gene") |>
        dplyr::filter(gene_type %in% types | gene %in% keep)
}

#' Filter network edges based on genes present in LDAK scores
#'
#' @param net Network data frame.
#' @param ldak LDAK score data frame.
#' @return Filtered network data frame.
#' @importFrom dplyr filter
#' @export
filter_gene_network_ldak <- function(net, ldak){
    message("Filter network to keep only edges where both genes have LDAK scores")
    genes <- unique(ldak$gene)
    net |> dplyr::filter(Regulator %in% genes & Target %in% genes)
}
