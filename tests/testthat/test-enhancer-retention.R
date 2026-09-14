test_that("annotation retention is the union of source-specific policies", {
  config <- list(enhancer_classification = list(
    promoter_window = "strand_aware_tss", class1_label = "enh_class1"))
  classes <- tibble::tibble(name = letters[1:5],
    enhancer_class = c("enh_class1", "enh_class2", "enh_class3", "enh_class2", "enh_class3"),
    interaction_type = c("HiC_PO", "JEME", "JEME", "HiC_PO", "HiC_PO;JEME"))
  expect_equal(retained_enhancer_annotations(classes, config)$name, c("a", "b", "c", "e"))
  config$enhancer_classification$promoter_window <- "gene_body"
  expect_equal(retained_enhancer_annotations(classes, config)$name, "a")
  config$enhancer_classification$drop_non_class1_edge_sources <- list()
  expect_equal(retained_enhancer_annotations(classes, config)$name, letters[1:5])
})
