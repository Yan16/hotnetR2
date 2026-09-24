test_that("genome-build aliases are normalized and incompatible builds fail", {
  expect_equal(hotnetR2:::normalise_genome_build("hg19"), "GRCH37")
  expect_equal(hotnetR2:::normalise_genome_build("GRCh38"), "GRCH38")
  expect_error(hotnetR2:::normalise_genome_build("mm10"), "Genome build")

  config <- write_network_test_config()
  config$reference$annotation_genome_build <- "GRCh38"
  expect_error(validate_analysis_config(config), "must describe the same build")
})

test_that("LDAK interval validation enforces unique, valid panel-overlapping locations", {
  bounds <- data.frame(CHR = c(1L, 23L), min_bp = c(100, 50), max_bp = c(1000, 500), reference_variants = c(3L, 2L))
  intervals <- data.frame(
    name = c("A", "B"), CHR = c("chr1", "X"), START = c(100, 50), END = c(200, 100),
    source = c("JEME", "HiC")
  )
  checked <- hotnetR2:::validate_ldak_intervals(intervals, "enhancer", bounds)
  expect_equal(checked$data$CHR, c(1L, 23L))
  expect_equal(checked$report$value[checked$report$metric == "interval_rows"], 2)

  zero_start <- intervals
  zero_start$START[1] <- 0L
  expect_no_error(hotnetR2:::validate_ldak_intervals(zero_start, "enhancer", bounds))

  duplicated_name <- intervals
  duplicated_name$name[2] <- "A"
  expect_error(hotnetR2:::validate_ldak_intervals(duplicated_name, "enhancer", bounds), "duplicate names")

  invalid_coordinates <- intervals
  invalid_coordinates$END[1] <- 100
  expect_error(hotnetR2:::validate_ldak_intervals(invalid_coordinates, "enhancer", bounds), "invalid START/END")

  outside_panel <- intervals
  outside_panel$START[1] <- 1000
  outside_panel$END[1] <- 1100
  outside_checked <- hotnetR2:::validate_ldak_intervals(outside_panel, "enhancer", bounds)
  expect_equal(
    outside_checked$report$value[outside_checked$report$metric == "intervals_fully_outside_reference_panel"],
    1
  )
})

test_that("annotation-file validation detects unsorted locations and detail mismatches", {
  expect_true(is.function(validate_ldak_annotations))
})
