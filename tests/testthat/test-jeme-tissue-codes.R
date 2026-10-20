test_that("JEME node provenance prefers stable tissue codes", {
  jeme <- tibble::tibble(
    file = c("lasso.90.bed.gz", "lasso.92.bed.gz"),
    nfile = c("90", "92"),
    tiss = c("E092", "E094"),
    desc1 = c("GI_STMC_FET__Fetal Stomach", "GI_STMC_GAST__Gastric"),
    enhancer = c("chr1:10-20", "chr1:10-20"),
    CHR = c("1", "1"),
    START = c(10L, 10L),
    END = c(20L, 20L),
    promoterFull = c("ENSG000001$P1", "ENSG000001$P1"),
    ENSG = c("ENSG000001", "ENSG000001"),
    promoter = c("P1", "P1"),
    CHR2 = c("1", "1"),
    location = c(100L, 100L)
  )

  nodes <- extract_jeme_nodes(jeme)

  expect_equal(nodes$enhancers$tissues, "E092, E094")
  expect_equal(nodes$promoters$tissues, "E092, E094")
  expect_false(grepl("Fetal Stomach|Gastric", nodes$enhancers$tissues))
})

test_that("JEME node provenance recovers codes from simplified-cache nfile", {
  jeme <- tibble::tibble(
    file = "lasso.90.bed.gz",
    nfile = "90",
    desc1 = "GI_STMC_FET__Fetal Stomach",
    enhancer = "chr1:10-20",
    CHR = "1",
    START = 10L,
    END = 20L,
    promoterFull = "ENSG000001$P1",
    ENSG = "ENSG000001",
    promoter = "P1",
    CHR2 = "1",
    location = 100L
  )

  nodes <- extract_jeme_nodes(jeme)

  expect_equal(nodes$enhancers$tissues, "E092")
  expect_equal(nodes$promoters$tissues, "E092")
})

test_that("JEME node provenance falls back for unknown legacy identifiers", {
  jeme <- tibble::tibble(
    file = "legacy.bed.gz",
    nfile = "9999",
    desc1 = "Legacy tissue description",
    enhancer = "chr1:10-20",
    CHR = "1",
    START = 10L,
    END = 20L,
    promoterFull = "ENSG000001$P1",
    ENSG = "ENSG000001",
    promoter = "P1",
    CHR2 = "1",
    location = 100L
  )

  nodes <- extract_jeme_nodes(jeme)

  expect_equal(nodes$enhancers$tissues, "Legacy tissue description")
  expect_equal(nodes$promoters$tissues, "Legacy tissue description")
})
