test_that("vignette analysis examples validate as current canonical profiles", {
  examples <- system.file(
    "extdata", "vignette_examples", package = "hotnetR2", mustWork = TRUE
  )
  files <- file.path(examples, c("analysis_aracne.yaml", "analysis_stringdb.yaml"))
  expect_true(all(file.exists(files)))

  configs <- lapply(files, read_analysis_config)
  expect_identical(
    vapply(configs, interaction_network_type, character(1)),
    c("aracne", "stringdb")
  )
  for (config in configs) {
    expect_false(config$enhancer_classification$enabled)
    expect_identical(config$ldak$promoter_interval_mode, "strand_aware_tss")
    expect_true(config$regulatory$promoter_harmonization$enabled)
    expect_identical(
      config$regulatory$promoter_harmonization$jeme_strategy,
      "all_tissues_ensg_map"
    )
    expect_identical(
      config$regulatory$promoter_harmonization$hic_strategy,
      "gene_symbol_alias"
    )
    expect_setequal(
      config$networks$selected,
      c("network_1", "network_2", "network_3", "network_4")
    )
  }
})

test_that("vignette launchers parse", {
  examples <- system.file(
    "extdata", "vignette_examples", package = "hotnetR2", mustWork = TRUE
  )
  expect_silent(parse(file.path(examples, "run_analysis.R")))
  testthat::skip_on_os("windows")
  status <- system2("bash", c("-n", file.path(examples, "run_analysis.sh")))
  expect_identical(status, 0L)
})
