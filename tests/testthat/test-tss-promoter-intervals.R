test_that("strand-aware TSS conversion uses start for plus and end for minus", {
  input <- data.frame(
    name = c("PLUS", "MINUS"), CHR = c(1L, 2L),
    START = c(100L, 300L), END = c(200L, 500L), strand = c("+", "-")
  )
  result <- hotnetR2:::convert_promoter_details_to_tss(input)
  expect_equal(result$TSS, c(101L, 500L))
  expect_equal(result$START, result$TSS - 1L)
  expect_equal(result$END, result$TSS)
  expect_equal(result$END - result$START, c(1L, 1L))
  expect_equal(result$gene_START, input$START)
  expect_equal(result$gene_END, input$END)
  expect_true(all(result$promoter_interval_mode == "strand_aware_tss"))
})

test_that("strand-aware TSS conversion rejects unknown strands", {
  input <- data.frame(name = "BAD", CHR = 1L, START = 1L, END = 2L, strand = "?")
  expect_error(
    hotnetR2:::convert_promoter_details_to_tss(input),
    "unsupported strand"
  )
})
