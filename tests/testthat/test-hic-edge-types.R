test_that("HiC edge types are validated and selected before cache access", {
  config <- list(regulatory = list(hic = list(tissue_type = "Gastric", edge_types = "PP")))
  local_mocked_bindings(get_hic_by_tissue = function(tissue_type, edge_type, cache_dir) {
    expect_equal(edge_type, "PP")
    expect_true(is.null(tissue_type) || identical(tissue_type, "Gastric"))
    tibble::tibble(Promoter1 = "A", Promoter2 = "B", Tissue_type = "Gastric")
  })
  expect_equal(nrow(load_analysis_hic(config)$PO), 0L)
  expect_equal(nrow(load_analysis_hic(config, all_tissues = TRUE)$PP), 1L)
  for (bad in list(list(), "P-P", c("PP", "PP"), NA_character_)) {
    config$regulatory$hic$edge_types <- bad
    expect_error(analysis_hic_edge_types(config), "edge_types")
  }
  config$regulatory$hic$edge_types <- NULL
  expect_equal(analysis_hic_edge_types(config), c("PO", "PP"))
  config$regulatory$hic$edge_types <- "PO"
  local_mocked_bindings(get_hic_by_tissue = function(tissue_type, edge_type, cache_dir) {
    expect_equal(edge_type, "PO")
    empty_analysis_hic()$PO
  })
  expect_equal(nrow(load_analysis_hic(config)$PP), 0L)
})
