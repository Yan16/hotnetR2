write_ldak_summary_fixture <- function(config, node_type) {
  paths <- analysis_paths(config)
  result_dir <- file.path(paths[["ldak_results"]], node_type)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  details_file <- file.path(paths[["ldak_annotations"]], paste0(node_type, "s.details.tsv.gz"))
  dir.create(dirname(details_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(data.frame(name = c("A", "B"), source = "fixture"), details_file)
  writeLines(c(
    "Gene_Name Heritability SE LRT_P_Perm",
    "A 0.1 0.2 0.01",
    "B 0.2 0.3 0.20"
  ), file.path(result_dir, "remls.all"))
  writeLines(c(
    "# LDAK fixture", "# generated for test", "Gene_Name Min_Pvalue",
    "A 0.001", "B 0.05"
  ), file.path(result_dir, "genes.details"))
}

test_that("summaries regenerate from completed LDAK files without rerunning LDAK", {
  config <- write_network_test_config()
  write_ldak_summary_fixture(config, "enhancer")
  write_ldak_summary_fixture(config, "promoter")

  first <- summarize_analysis_ldak_results(config, dry_run = FALSE)
  enhancer_file <- ldak_summary_specs(config)$enhancer$output_file
  first_table <- readr::read_tsv(enhancer_file, show_col_types = FALSE)
  second <- summarize_analysis_ldak_results(config, dry_run = FALSE)
  second_table <- readr::read_tsv(enhancer_file, show_col_types = FALSE)

  expect_equal(first$reml_rows, c(2, 2))
  expect_equal(second$unmatched_annotations, c(0, 0))
  expect_equal(first_table, second_table)
  expect_true(all(c("SE", "SD", "Min_Pvalue", "FDR", "gene", "flank") %in% names(first_table)))
  expect_equal(first_table$SD, first_table$SE)
  expect_equal(first_table$Min_Pvalue, c(0.001, 0.05))
})

test_that("standardized-score comparison quantifies common score fields", {
  observed_file <- tempfile(fileext = ".tsv")
  reference_file <- tempfile(fileext = ".tsv")
  readr::write_tsv(data.frame(Gene_Name = c("A", "B"), LRT_P_Perm = c(0.01, 0.20), FDR = c(0.02, 0.20)), observed_file)
  readr::write_tsv(data.frame(Gene_Name = c("A", "B"), LRT_P_Perm = c(0.01, 0.25), FDR = c(0.02, 0.25)), reference_file)
  comparison <- hotnetR2:::compare_ldak_summary_tables(observed_file, reference_file)
  pvalues <- comparison[comparison$metric == "LRT_P_Perm", ]
  expect_equal(pvalues$shared_ids, 2)
  expect_equal(pvalues$different_values, 1)
  expect_equal(pvalues$max_abs_difference, 0.05)
})
