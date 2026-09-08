test_that("provenance files and event logs remain below the analysis directory", {
  config <- write_network_test_config()
  files <- write_analysis_provenance(config, checksum = FALSE)
  root_prefix <- paste0(normalizePath(dirname(config$config_file)), .Platform$file.sep)

  expect_true(all(file.exists(files)))
  expect_true(all(startsWith(unname(files), root_prefix)))
  manifest <- readr::read_tsv(files[["input_manifest"]], show_col_types = FALSE)
  expect_true(all(is.na(manifest$md5)))
  events <- readr::read_tsv(files[["event_log"]], show_col_types = FALSE)
  expect_true(any(events$event == "provenance_written"))
})
