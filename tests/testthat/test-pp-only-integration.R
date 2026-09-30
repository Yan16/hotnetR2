test_that("PP-only annotations and networks exclude PO while retaining JEME and PP", {
  skip_if(Sys.which("bedtools") == "")
  config <- write_network_test_config(available = "[network_1, network_2]", selected = "[network_1, network_2]")
  config$aracne <- list(enabled = FALSE)
  config$regulatory$jeme <- list(key_column = "tiss", value = "E094")
  config$regulatory$hic <- list(tissue_type = "Gastric", edge_types = "PP")
  config$regulatory$promoter_harmonization <- list(enabled = TRUE, gtf_file = "gtf.rds", hic_alias_file = "aliases.rds")
  config$reference$ncbi_feature_table <- "ncbi.txt"
  config$networks$split_multi_promoter_edges <- TRUE
  config$networks$deduplicate_jeme_edges_before_harmonization <- TRUE
  config$enhancer_classification <- list(enabled = TRUE, source_analysis = ".", enhancer_flank_bp = 0L,
    promoter_flank_bp = 0L, promoter_window = "strand_aware_tss",
    class1_label = "enh_class1", class2_label = "enh_class2", class3_label = "enh_class3")
  config$ldak$promoter_interval_mode <- "strand_aware_tss"
  config$ldak$annotation_tissues <- list(jeme = "ALL", hic = "ALL")
  readr::write_lines("fixture", file.path(config$project_root, "ncbi.txt"))
  saveRDS(tibble::tibble(gene_id = "ENSG1", gene_name = "A"), file.path(config$project_root, "gtf.rds"))
  saveRDS(tibble::tibble(gene = c("A", "B"), Hsym = c("A", "B")), file.path(config$project_root, "aliases.rds"))
  jeme <- tibble::tibble(enhancer = "chr1:100-110", CHR = "1", START = 100L, END = 110L,
    promoterFull = "A", promoter = "A", ENSG = "ENSG1", CHR2 = "1", location = 500L,
    desc1 = "Gastric", nfile = 92L)
  local_mocked_bindings(
    get_jeme = function(...) jeme,
    get_hic = function(...) stop("Both-source reader must not be called"),
    get_hic_by_tissue = function(tissue_type, edge_type, cache_dir) {
      if (edge_type != "PP") stop("PO cache must not be read")
      tibble::tibble(Promoter1 = "A", Promoter2 = "B", Tissue_type = "Gastric")
    },
    NCBI_feature_table = function(...) tibble::tibble(gene_symbol = c("A", "B"), gene_name = c("A", "B"),
      CHR = 1L, START = c(500L, 700L), END = c(550L, 750L), strand = "+"),
    ldak_reference_bounds = function(...) tibble::tibble(CHR = 1L, min_bp = 1, max_bp = 1000, reference_variants = 2L)
  )
  expect_no_error(validate_analysis_config(config))
  setup_analysis(config, dry_run = FALSE)
  expect_no_error(prepare_ldak_annotations(config, dry_run = FALSE))
  classes <- classify_enhancers(config)
  expect_equal(classes$source, "JEME")
  expect_equal(classes$enhancer_class, "enh_class1")
  paths <- analysis_paths(config)
  details <- readr::read_tsv(file.path(paths[["ldak_annotations"]], "promoters.details.tsv.gz"), show_col_types = FALSE)
  expect_setequal(details$name, c("A", "B"))
  expect_false(any(grepl("HiC_PO", details$interaction_type)))
  for (kind in c("enhancer", "promoter")) {
    readr::write_tsv(tibble::tibble(Gene_Name = if (kind == "enhancer") "chr1:100-110" else c("A", "B"),
      cohort = "test", LRT_P_Perm = 1e-6), file.path(paths[["ldak_summary"]], paste0(kind, "_ldak.tsv.gz")))
  }
  expect_no_error(build_hhotnet_networks(config, dry_run = FALSE))
  edges <- readr::read_tsv(file.path(paths[["hhotnet_data"]], "regulatory_edges.tsv"), show_col_types = FALSE)
  expect_setequal(edges$source, c("JEME", "HiC_PP"))
  expect_equal(nrow(edges), 2L)
})
