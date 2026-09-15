test_that("curated aliases resolve both conflicts and preserve unique lookup keys", {
  tables <- packaged_alias_tables()
  expect_equal(nrow(tables$alias), 60422L)
  expect_equal(nrow(tables$alias_nodup), 59241L)
  expect_false(anyDuplicated(tables$alias_nodup$gene) > 0L)
  selected <- dplyr::filter(tables$alias_nodup, .data$gene %in% c("C11orf48", "MEGT1"))
  expect_equal(selected$Hsym, c("LBHD1", "LY6G6D"))
  expect_equal(selected$Hhgnc_id, c("HGNC:28351", "HGNC:13935"))
  expect_equal(attr(tables$alias_nodup, "mapping_direction"), "raw_to_approved")
})

test_that("cache generation uses packaged references and preserves existing files", {
  cache <- withr::local_tempdir()
  result <- create_alias_references(cache)
  expect_true(all(result$status == "generated"))
  expected <- packaged_alias_tables()
  expect_identical(readRDS(result$path[[1]]), expected$alias)
  expect_identical(readRDS(result$path[[2]]), expected$alias_nodup)
  saveRDS("sentinel", result$path[[1]])
  expect_true(all(create_alias_references(cache)$status == "skipped_existing"))
  expect_identical(readRDS(result$path[[1]]), "sentinel")
  create_alias_references(cache, overwrite = TRUE)
  expect_identical(readRDS(result$path[[1]]), expected$alias)
})

test_that("full reference setup uses the same curated default without rewriting cached GENCODE", {
  skip_if_not_installed("rtracklayer")
  cache <- withr::local_tempdir()
  dir.create(file.path(cache, "GTEx"))
  readr::write_lines("not imported because GENCODE is already cached",
    file.path(cache, "GTEx", "gencode.v39.GRCh38.genes.gtf"))
  gencode <- file.path(cache, "gencode.v39.GRCh38.genes_nodup.rds")
  saveRDS(tibble::tibble(gene_name = "sentinel"), gencode)
  result <- create_harmonization_references(cache)
  expect_equal(result$status, c("skipped_existing", "generated", "generated"))
  expect_identical(readRDS(result$path[[3]]), packaged_alias_tables()$alias_nodup)
  expect_equal(readRDS(gencode)$gene_name, "sentinel")
})

test_that("HiC maps raw labels forward with new files but keeps legacy direction", {
  file <- withr::local_tempfile(fileext = ".rds")
  aliases <- tibble::tibble(gene = c("C11orf48", "MEGT1"), Hsym = c("LBHD1", "LY6G6D"))
  attr(aliases, "mapping_direction") <- "raw_to_approved"
  saveRDS(aliases, file)
  hic <- list(PO = tibble::tibble(Regulator = "enh", Target = "C11orf48"),
              PP = tibble::tibble(Promoter1 = "MEGT1", Promoter2 = "C11orf48"))
  mapped <- harmonize_hic_gtf(hic, file, verbose = FALSE)
  expect_equal(mapped$PO$Target_harmonized, "LBHD1")
  expect_equal(mapped$PP$Promoter1_harm, "LY6G6D")
  expect_equal(mapped$PP$Promoter2_harm, "LBHD1")
  aliases$ensg <- c("ENSG0001", "ENSG0002")
  saveRDS(aliases, file)
  candidates <- aracne_mapping_candidates(list(map_ensembl = FALSE,
    map_aliases = TRUE, alias_file = file), dirname(file))
  mapped_aracne <- map_aracne_identifiers(c("C11orf48", "MEGT1", "ENSG0001"), candidates)
  expect_equal(mapped_aracne$canonical, c("LBHD1", "LY6G6D", "LBHD1"))
  attr(aliases, "mapping_direction") <- NULL
  saveRDS(aliases, file)
  hic$PO$Target <- "LBHD1"
  expect_equal(harmonize_hic_gtf(hic, file, verbose = FALSE)$PO$Target_harmonized, "C11orf48")
})

test_that("external builder reproduces packaged data from unchanged workbooks", {
  builder <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "harmonization", "build_alias_data.R", package = "hotnetR2"), envir = builder)
  output <- withr::local_tempdir()
  rebuilt <- builder$build_alias_data(
    system.file("extdata", "harmonization_2026_0909", package = "hotnetR2"),
    system.file("extdata", "gtex_hugo_merged_by_ens_SORTens_FINALv1.xlsx", package = "hotnetR2"), output)
  packaged <- packaged_alias_tables()
  expect_identical(rebuilt$alias, packaged$alias)
  expect_identical(rebuilt$alias_nodup, packaged$alias_nodup)
  audit <- readr::read_tsv(file.path(output, "manual_resolution_audit.tsv"), show_col_types = FALSE)
  expect_setequal(audit$gene, c("C11orf48", "MEGT1"))
  expect_true(all(audit$Hsym[audit$gene == "C11orf48"] == "LBHD1"))
  expect_true(all(audit$Hsym[audit$gene == "MEGT1"] == "LY6G6D"))
})
