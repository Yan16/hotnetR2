test_that("local Python h-HotNet runner uses the configured source and compatibility shim", {
  config <- write_network_test_config(
    available = "[network_1, network_2, network_3, network_4]",
    selected = "[network_1, network_2, network_3, network_4]"
  )
  config$hhotnet$execution_mode <- "local_python"
  config$hhotnet$local_dir <- "/tmp/hierarchical-hotnet"
  config$hhotnet$python_executable <- "python3"
  config$hhotnet$compat_dir <- "compat"
  config$hhotnet$apptainer_image <- "unused.sif"
  lines <- make_hhotnet_job_lines(config, networks = "network_4")
  expect_true(any(grepl("DEFAULT_MODE='local_python'", lines, fixed = TRUE)))
  expect_true(any(grepl("${HHOTNET_DIR}/src/construct_hierarchy.py", lines, fixed = TRUE)))
  expect_true(any(grepl("DEFAULT_NETWORKS=(network_4)", lines, fixed = TRUE)))
  expect_true(any(grepl("Network is already running", lines, fixed = TRUE)))
  script <- tempfile(fileext = ".sh")
  writeLines(lines, script)
  expect_equal(system2("bash", c("-n", script)), 0)
})
