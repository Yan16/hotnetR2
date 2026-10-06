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
  local_mocked_bindings(
    get_jeme = function(...) tibble::tibble(
      enhancer = c("A", "A", "A"), promoter = c("P1", "P2", "P1"),
      promoterFull = c("ENSG000001$P1", "ENSG000002$P2", "ENSG000001$P1"),
      ENSG = c("ENSG000001", "ENSG000002", "ENSG000001"),
      conf_score = c(0.91, 0.72, 0.91),
      nfile = c("90", "90", "90")
    )
  )

  first <- summarize_analysis_ldak_results(config, dry_run = FALSE)
  enhancer_file <- ldak_summary_specs(config)$enhancer$output_file
  long_file <- enhancer_ldak_long_file(config)
  first_table <- readr::read_tsv(enhancer_file, show_col_types = FALSE)
  promoter_table <- readr::read_tsv(
    ldak_summary_specs(config)$promoter$output_file, show_col_types = FALSE
  )
  long_table <- readr::read_tsv(long_file, show_col_types = FALSE)
  second <- summarize_analysis_ldak_results(config, dry_run = FALSE)
  second_table <- readr::read_tsv(enhancer_file, show_col_types = FALSE)

  expect_equal(first$reml_rows, c(2, 2))
  expect_equal(second$unmatched_annotations, c(0, 0))
  expect_equal(first_table, second_table)
  expect_true(all(c("Gene_Name", "original_promoter", "SE", "SD",
                    "Min_Pvalue", "FDR", "flank") %in% names(first_table)))
  expect_false("gene" %in% names(first_table))
  expect_false("gene" %in% names(promoter_table))
  expect_equal(first_table$original_promoter, c("P1;P2", NA_character_))
  expect_equal(first_table$SD, first_table$SE)
  expect_equal(first_table$Min_Pvalue, c(0.001, 0.05))
  expect_equal(nrow(long_table), 3L)
  expect_equal(long_table$promoterFull, c("ENSG000001$P1", "ENSG000002$P2", NA_character_))
  expect_equal(long_table$ENSG, c("ENSG000001", "ENSG000002", NA_character_))
  expect_equal(long_table$original_promoter, c("P1", "P2", NA_character_))
  expect_equal(long_table$promoter, c("P1", "P2", NA_character_))
  expect_equal(long_table$conf_score, c(0.91, 0.72, NA_real_))
  expect_equal(long_table$tissue, c("E092", "E092", NA_character_))
  expect_equal(long_table$tissue_name, c("Fetal Stomach", "Fetal Stomach", NA_character_))
})

test_that("enhancer long summary expands unique targets across JEME tissues", {
  config <- write_network_test_config()
  output_dir <- analysis_paths(config)[["ldak_summary"]]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(
    tibble::tibble(Gene_Name = c("E1", "E2"), gene = c("E1", "E2"),
                   LRT_P_Perm = c(0.01, 0.2)),
    ldak_summary_specs(config)$enhancer$output_file
  )
  local_mocked_bindings(
    get_jeme = function(...) tibble::tibble(
      enhancer = c("E1", "E1", "E1", "E1", "E3"),
      promoter = c("P1", "P2", "P1", "P1", "P3"),
      promoterFull = c("ENSG000001$P1", "ENSG000002$P2", "ENSG000001$P1", "ENSG000001$P1", "ENSG000003$P3"),
      ENSG = c("ENSG000001", "ENSG000002", "ENSG000001", "ENSG000001", "ENSG000003"),
      conf_score = c(0.91, 0.72, 0.88, 0.88, 0.65),
      nfile = c("90", "90", "92", "92", "90")
    )
  )

  output <- create_enhancer_ldak_long_summary(config)
  observed <- readr::read_tsv(output, show_col_types = FALSE)
  standard <- readr::read_tsv(
    ldak_summary_specs(config)$enhancer$output_file, show_col_types = FALSE
  )

  expect_equal(basename(output), "enhancer_ldak_long.tsv.gz")
  expect_equal(nrow(observed), 4L)
  expect_equal(observed$Gene_Name, c("E1", "E1", "E1", "E2"))
  expect_equal(observed$promoterFull,
               c("ENSG000001$P1", "ENSG000002$P2", "ENSG000001$P1", NA_character_))
  expect_equal(observed$ENSG,
               c("ENSG000001", "ENSG000002", "ENSG000001", NA_character_))
  expect_equal(observed$original_promoter, c("P1", "P2", "P1", NA_character_))
  expect_equal(observed$promoter, c("P1", "P2", "P1", NA_character_))
  expect_equal(observed$conf_score, c(0.91, 0.72, 0.88, NA_real_))
  expect_equal(observed$tissue, c("E092", "E092", "E094", NA_character_))
  expect_equal(observed$tissue_name,
               c("Fetal Stomach", "Fetal Stomach", "Gastric", NA_character_))
  expect_false("gene" %in% names(observed))
  expect_false("gene" %in% names(standard))
  expect_equal(standard$original_promoter, c("P1;P2", NA_character_))
})

test_that("enhancer summary rejects a conflicting legacy gene alias", {
  config <- write_network_test_config()
  output_dir <- analysis_paths(config)[["ldak_summary"]]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(
    tibble::tibble(Gene_Name = "E1", gene = "different", LRT_P_Perm = 0.01),
    ldak_summary_specs(config)$enhancer$output_file
  )
  local_mocked_bindings(
    get_jeme = function(...) tibble::tibble(
      enhancer = "E1", promoter = "P1", promoterFull = "ENSG000001$P1",
      ENSG = "ENSG000001", conf_score = 0.91, nfile = "90"
    )
  )

  expect_error(
    create_enhancer_ldak_long_summary(config),
    "Legacy gene column disagrees with Gene_Name"
  )
})

test_that("enhancer long summary rejects conflicting association scores", {
  config <- write_network_test_config()
  output_dir <- analysis_paths(config)[["ldak_summary"]]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(
    tibble::tibble(Gene_Name = "E1", LRT_P_Perm = 0.01),
    ldak_summary_specs(config)$enhancer$output_file
  )
  local_mocked_bindings(
    get_jeme = function(...) tibble::tibble(
      enhancer = c("E1", "E1"), promoter = c("P1", "P1"),
      promoterFull = c("ENSG000001$P1", "ENSG000001$P1"),
      ENSG = c("ENSG000001", "ENSG000001"),
      conf_score = c(0.91, 0.72), nfile = c("90", "90")
    )
  )

  expect_error(
    create_enhancer_ldak_long_summary(config),
    "conflicting conf_score values"
  )
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
