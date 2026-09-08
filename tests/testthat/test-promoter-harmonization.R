test_that("all-tissue JEME map harmonizes by versionless ENSG", {
  all_jeme <- data.frame(
    ENSG = c("ENSG000001.1", "ENSG000001.2", "ENSG000002.1", "ENSG000003.1"),
    promoter = c("OLD1", "OLDER1", "GENE2", "UNCHANGED")
  )
  gtf <- data.frame(
    gene_id = c("ENSG000001.9", "ENSG000002.4"),
    gene_name = c("GENE1", "GENE2")
  )
  mapping <- hotnetR2:::build_all_tissue_jeme_ensg_map(all_jeme, gtf)
  expect_equal(nrow(mapping), 3L)
  expect_equal(mapping$promoter_harmonized[mapping$ensg_key == "ENSG000001"], "GENE1")
  expect_equal(mapping$input_promoter_count[mapping$ensg_key == "ENSG000001"], 2L)

  selected <- data.frame(
    ENSG = c("ENSG000001.7", "ENSG000003.2", "ENSG999999.1"),
    promoter = c("OLD1", "UNCHANGED", "NOVEL"),
    enhancer = c("e1", "e2", "e3")
  )
  result <- hotnetR2:::apply_all_tissue_jeme_ensg_map(selected, mapping)
  expect_equal(result$promoter, c("GENE1", "UNCHANGED", "NOVEL"))
  expect_equal(result$enhancer, selected$enhancer)
  expect_equal(
    result$harmonization_status,
    c(
      "mapped_by_ensg", "unmapped_ensg_fallback_to_input",
      "ensg_absent_from_all_tissue_map_fallback_to_input"
    )
  )
})

test_that("ambiguous GENCODE ENSG mappings fall back to supplied JEME symbols", {
  all_jeme <- data.frame(ENSG = "ENSG000001.1", promoter = "INPUT")
  gtf <- data.frame(
    ENSG = c("ENSG000001", "ENSG000001"), gene = c("GENEA", "GENEB")
  )
  mapping <- hotnetR2:::build_all_tissue_jeme_ensg_map(all_jeme, gtf)
  result <- hotnetR2:::apply_all_tissue_jeme_ensg_map(all_jeme, mapping)
  expect_equal(result$promoter, "INPUT")
  expect_equal(result$harmonization_status, "ambiguous_ensg_fallback_to_input")
})
