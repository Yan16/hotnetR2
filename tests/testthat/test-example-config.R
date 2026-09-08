test_that("the shipped example configuration validates", {
  example <- system.file(
    "extdata", "analysis.example.yaml",
    package = "hotnetR2", mustWork = TRUE
  )
  config <- read_analysis_config(example)
  expect_equal(config$schema_version, 1)
  expect_equal(available_networks(config), paste0("network_", 1:4))
  expect_equal(resolve_networks(config), paste0("network_", 1:3))
  expect_equal(
    unname(analysis_paths(config)[["root"]]),
    normalizePath(dirname(example), mustWork = TRUE)
  )
})
