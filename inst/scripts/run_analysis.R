#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript run_analysis.R analysis.yaml [stage ...]")
if (length(args) == 1L) {
  hotnetR2::run_analysis(args[[1]])
} else {
  hotnetR2::run_analysis(args[[1]], stages = args[-1L])
}
