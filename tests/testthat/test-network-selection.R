test_that("network resolution supports defaults, a single network, and subsets", {
  config <- write_network_test_config()
  expect_equal(resolve_networks(config), c("network_1", "network_3"))
  expect_equal(resolve_networks(config, "network_2"), "network_2")
  expect_equal(resolve_networks(config, c("network_3", "network_1")), c("network_3", "network_1"))
  expect_error(resolve_networks(config, "network_4"), "not available")
  expect_error(resolve_networks(config, c("network_1", "network_1")), "duplicate")
})

test_that("selected networks must be a subset of available networks", {
  expect_error(
    write_network_test_config(available = "[network_1]", selected = "[network_2]"),
    "must occur in networks.available"
  )
})

test_that("configured input paths are reported without needing an external file", {
  config <- write_network_test_config()
  report <- check_analysis_inputs(config, write_report = FALSE)
  expect_true(all(report$status == "FAIL"))
  expect_true(all(report$detail == "missing"))
})

test_that("legacy and current edge-filter modes are accepted but unknown modes fail", {
  config <- write_network_test_config()
  config$networks$regulatory_edge_filter <- "legacy_filter_network_ldak"
  expect_silent(validate_analysis_config(config))
  config$networks$regulatory_edge_filter <- "unsupported_mode"
  expect_error(validate_analysis_config(config), "regulatory_edge_filter")
})

test_that("effective configuration is analysis-local and records resolved networks", {
  config <- write_network_test_config()
  file <- write_effective_config(config)
  saved <- yaml::read_yaml(file)
  expect_true(file.exists(file))
  expect_equal(unlist(saved$networks$available), c("network_1", "network_2", "network_3"))
  expect_equal(unlist(saved$networks$selected), c("network_1", "network_3"))
  expect_error(write_effective_config(config, tempfile("outside-")), "below the analysis directory")
})
