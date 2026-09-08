#' hotnetR: Analysis of Biological Networks with HotNet
#'
#' @description
#' The `hotnetR` package provides a collection of functions for biological
#' network analysis, with a focus on data preparation for HotNet, GTEx, and LDAK studies.
#' It includes utilities for downloading interaction data (JEME, HiC),
#' harmonizing genomic references, and preparing SLURM-compatible scripts
#' for downstream tools like ARACNe.
#'
#' @section Key Functions:
#' \itemize{
#'   \item \code{\link{setup_data}}: Download all required genomic datasets to the local cache.
#'   \item \code{\link{get_jeme}}: Access JEME interaction datasets by tissue.
#'   \item \code{\link{get_hic}}: Access processed HiC map datasets.
#'   \item \code{\link{export_gtex_aracne}}: Prepare GTEx expression data for ARACNe analysis.
#'   \item \code{\link{create_aracne_slurm}}: Generate SLURM scripts for parallelized network construction.
#' }
#'
"_PACKAGE"
