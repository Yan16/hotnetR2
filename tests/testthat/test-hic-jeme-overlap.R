test_that("HiC-JEME overlaps use both flanks, chromosomes and inclusive boundaries", {
  po <- tibble::tibble(Interacting_fragment = c("hit", "adjacent", "other_chr", "sex", "hit"),
                       CHR = c("1", "1", "2", "chrX", "1"),
                       START = c(3000, 3001, 3000, 3000, 3000),
                       END = c(3010, 3010, 3010, 3010, 3010))
  jeme <- tibble::tibble(enhancer = c("j1", "jx"), CHR = c("chr1", "23"), START = 900, END = 1000)
  result <- hic_jeme_overlap_tables(po, jeme)
  expect_setequal(result$overlaps$hic_enhancer, c("hit", "sex"))
  expect_equal(nrow(result$classification), 4L)
  expect_setequal(dplyr::filter(result$classification, .data$retained_by_jeme_filter)$hic_enhancer,
                  c("adjacent", "other_chr"))
  expect_equal(result$overlaps$hic_start, c(2000, 2000))
  expect_equal(nrow(hic_jeme_overlap_tables(po, jeme, 0, 0)$overlaps), 0L)
  expect_true(all(hic_jeme_overlap_tables(po, jeme[0, ])$classification$retained_by_jeme_filter))
  expect_equal(nrow(hic_jeme_overlap_tables(po[0, ], jeme)$classification), 0L)
})

test_that("both source stages use selected tissues and preserve PP and unfiltered defaults", {
  root <- withr::local_tempdir()
  config <- list(cache_dir = root, regulatory = list(
    jeme = list(key_column = "tiss", value = c("E094", "E059")),
    hic = list(tissue_type = "Gastric", jeme_overlap = list(enabled = TRUE,
      hic_flank_bp = 1000, jeme_flank_bp = 1000))))
  po <- tibble::tibble(Interacting_fragment = c("hit", "keep", "hit"),
    Promoter = c("A", "B", "C"), CHR = "1", START = c(3000, 9000, 3000),
    END = c(3010, 9010, 3010), Tissue_type = "Gastric")
  pp <- tibble::tibble(Promoter1 = "A", Promoter2 = "B", Tissue_type = "Gastric")
  local_mocked_bindings(
    get_hic = function(tissue_type, cache_dir) list(PO = po, PP = pp),
    get_jeme = function(key_column, value, method, simplified, cache_dir) {
      expect_identical(key_column, "tiss")
      expect_identical(value, c("E094", "E059"))
      # An unselected tissue's enhancer at 9000 must not enter this reference.
      list(E094 = tibble::tibble(enhancer = "j1", CHR = "1", START = 900, END = 1000),
           E059 = tibble::tibble(enhancer = "j2", CHR = "2", START = 9000, END = 9010))
    },
    annotation_output_file = function(config, name) file.path(root, name)
  )
  for (all_tissues in c(TRUE, FALSE)) {
    result <- load_analysis_hic(config, all_tissues)
    expect_identical(result$PP, pp)
    expect_identical(result$PO, po[2, ])
  }
  expect_length(list.files(root, pattern = "tsv.gz$"), 6L)
  config$regulatory$hic$jeme_overlap$enabled <- FALSE
  expect_identical(load_analysis_hic(config)$PO, po)
  config$regulatory$hic$jeme_overlap$hic_flank_bp <- -1
  expect_error(validate_hic_jeme_overlap(config), "nonnegative integer")
})
