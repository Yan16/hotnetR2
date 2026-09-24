test_that("strand-aware TSS conversion uses start for plus and end for minus", {
  input <- data.frame(
    name = c("PLUS", "MINUS"), CHR = c(1L, 2L),
    START = c(100L, 300L), END = c(200L, 500L), strand = c("+", "-")
  )
  result <- hotnetR2:::convert_promoter_details_to_tss(input)
  expect_equal(result$TSS, c(100L, 500L))
  expect_equal(result$START, result$TSS - 1L)
  expect_equal(result$END, result$TSS)
  expect_equal(result$END - result$START, c(1L, 1L))
  expect_equal(result$gene_START, input$START)
  expect_equal(result$gene_END, input$END)
  expect_true(all(result$promoter_interval_mode == "strand_aware_tss"))
})

test_that("NCBI gene bodies convert from 1-based inclusive to BED", {
  input <- data.frame(
    name = c("A", "B"), CHR = 1L,
    START = c(1L, 100L), END = c(10L, 200L), strand = c("+", "-")
  )
  result <- hotnetR2:::convert_promoter_details_to_gene_body(input)
  expect_equal(result$START, c(0L, 99L))
  expect_equal(result$END, c(10L, 200L))
  expect_equal(result$gene_START, input$START)
  expect_equal(result$gene_END, input$END)
  expect_true(all(result$promoter_interval_mode == "gene_body"))
})

test_that("promoter_range converts NCBI and JEME TSS values to one-base BED", {
  ncbi <- data.frame(
    chromosome = c("1", "2"), start = c(100L, 300L), end = c(200L, 500L),
    strand = c("+", "-")
  )
  ncbi_result <- promoter_range(ncbi, "NCBI")
  expect_equal(ncbi_result$START, c(99L, 499L))
  expect_equal(ncbi_result$END, c(100L, 500L))

  jeme <- data.frame(CHR2 = "chr1", location = 100L)
  jeme_result <- promoter_range(jeme, "jeme")
  expect_equal(jeme_result$START, 99L)
  expect_equal(jeme_result$END, 100L)
})

test_that("jeme_promoter_to_loc writes canonical one-base BED coordinates", {
  output_jeme <- withr::local_tempfile(fileext = ".loc")
  promoters <- data.frame(
    promoter = c("PLUS", "MINUS"), CHR = c(1L, 2L),
    promoterFull = c("p1", "p2"), ENSG = c("e1", "e2"),
    location = c(100L, 500L), strand = c("+", "-")
  )
  from_jeme <- jeme_promoter_to_loc(promoters, output_jeme, "jeme")
  expect_equal(from_jeme$START, c(99L, 499L))
  expect_equal(from_jeme$END, c(100L, 500L))

  output_ncbi <- withr::local_tempfile(fileext = ".loc")
  ncbi <- data.frame(
    gene_symbol = c("PLUS", "MINUS"), CHR = c(1L, 2L),
    START = c(100L, 300L), END = c(200L, 500L), strand = c("+", "-"),
    GeneID = c("1", "2")
  )
  from_ncbi <- jeme_promoter_to_loc(promoters, output_ncbi, "NCBI", ncbi)
  expect_equal(from_ncbi$START, c(99L, 499L))
  expect_equal(from_ncbi$END, c(100L, 500L))
})

test_that("standalone LDAK exporters preserve BED enhancers and convert NCBI TSS", {
  enhancer_file <- withr::local_tempfile(fileext = ".loc")
  nodes <- list(enhancers = data.frame(
    enhancer = "chr1:0-10", CHR = "chr1", START = 0L, END = 10L
  ))
  enhancer <- export_enhancer_ldak(nodes, enhancer_file, "jeme")
  expect_equal(enhancer$START, 0L)
  expect_equal(enhancer$END, 10L)

  root <- withr::local_tempdir()
  feature_file <- file.path(root, "feature_table.tsv")
  readr::write_tsv(data.frame(
    feature = c("gene", "gene"), chromosome = c("1", "2"),
    start = c(100L, 300L), end = c(200L, 500L),
    strand = c("+", "-"), symbol = c("PLUS", "MINUS"),
    check.names = FALSE
  ), feature_file)
  promoter_file <- file.path(root, "promoters.loc")
  create_ldak_promoter_boundary(feature_file, root, promoter_file)
  promoter <- readr::read_tsv(
    promoter_file,
    col_names = c("name", "CHR", "START", "END", "strand"),
    show_col_types = FALSE
  )
  expect_equal(promoter$START, c(99L, 499L))
  expect_equal(promoter$END, c(100L, 500L))
})

test_that("strand-aware TSS conversion rejects unknown strands", {
  input <- data.frame(name = "BAD", CHR = 1L, START = 1L, END = 2L, strand = "?")
  expect_error(
    hotnetR2:::convert_promoter_details_to_tss(input),
    "unsupported strand"
  )
})
