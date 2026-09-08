test_that("resource restoration is checksummed and independent of working directory", {
  source <- withr::local_tempdir(pattern = "bundle with spaces ")
  target <- withr::local_tempdir(pattern = "restored with spaces ")
  writeLines("fixture", file.path(source, "input.txt"))
  manifest <- capture_resources(source, "input.txt")
  provision_resources(manifest, target, source_root = source)
  expect_true(compare_artifacts(source, target, "input.txt")$identical)
  writeLines("changed", file.path(source, "input.txt"))
  expect_error(provision_resources(manifest, target, source_root = source), NA)
  fresh <- withr::local_tempdir()
  expect_error(provision_resources(manifest, fresh, source_root = source), "checksum")
  expect_false(file.exists(file.path(fresh, "input.txt")))
  manifest$path <- "../escape.txt"
  expect_error(provision_resources(manifest, fresh, source_root = source), "relative paths")
})

test_that("installed templates preserve both profiles without machine paths", {
  for (profile in c("v2", "v3")) {
    root <- withr::local_tempdir(pattern = "new project ")
    file <- initialize_analysis(root, profile)
    config <- read_analysis_config(file)
    expect_identical(config$enhancer_classification$promoter_window,
                     if (profile == "v2") "gene_body" else "strand_aware_tss")
    expect_identical(resolve_networks(config), c("network_1", "network_2"))
    expect_true(file.exists(file.path(root, "tools", "compat", "sitecustomize.py")))
    expect_false(any(grepl("/Users/", readLines(file), fixed = TRUE)))
    expect_error(initialize_analysis(root, profile), "already exists")
  }
})
