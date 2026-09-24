test_that("classification preserves containment precedence and chromosome separation", {
  skip_if(Sys.which("bedtools") == "", "bedtools not installed")
  config <- write_network_test_config()
  setup_analysis(config, dry_run = FALSE)
  config$enhancer_classification <- list(
    enhancer_flank_bp = 0L, promoter_flank_bp = 0L,
    promoter_window = "strand_aware_tss",
    class1_label = "enh_class1", class2_label = "enh_class2", class3_label = "enh_class3"
  )
  enhancer <- tibble::tibble(
    name = c("contained", "partial", "adjacent", "outside", "other_chr"),
    CHR = c(1L, 1L, 1L, 1L, 2L), START = c(12, 19, 20, 31, 12), END = c(14, 25, 21, 32, 14),
    source = "JEME", interaction_type = "JEME", tissues = "E094"
  )
  promoter <- tibble::tibble(name = c("A", "B"), CHR = 1L,
                              START = c(10, 13), END = c(20, 15))
  annotations <- analysis_paths(config)[["ldak_annotations"]]
  readr::write_tsv(enhancer, file.path(annotations, "enhancers.details.tsv.gz"))
  readr::write_tsv(promoter, file.path(annotations, "promoters.details.tsv.gz"))
  classes <- classify_enhancers(config)
  expect_identical(classes$enhancer_class,
                   c("enh_class2", "enh_class3", "enh_class1", "enh_class1", "enh_class1"))
  expect_identical(classes$fully_containing_promoter_count, c(1L, 0L, 0L, 0L, 0L))
  expect_identical(classes$overlapping_promoters[[1]], "A;B")
  expect_equal(nrow(readr::read_tsv(file.path(annotations, "enhancers.loc"),
                                  col_names = FALSE, show_col_types = FALSE)), 5)
  expect_equal(nrow(readr::read_tsv(file.path(annotations, "enhancers_class1.loc"),
                                  col_names = FALSE, show_col_types = FALSE)), 3)
  expect_identical(classify_enhancers(config), classes)
  config$enhancer_classification$drop_non_class1_edge_sources <- c("JEME", "HiC_PO")
  expect_identical(classify_enhancers(config), classes)
  expect_equal(nrow(readr::read_tsv(file.path(annotations, "enhancers.loc"),
                                  col_names = FALSE, show_col_types = FALSE)), 3)
})
