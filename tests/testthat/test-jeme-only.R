test_that("empty HiC selections disable loading, while named tissues retain behavior", {
  config <- write_network_test_config()
  for (selection in list(NULL, character(), list(), "", "  ", NA_character_)) {
    config$regulatory$hic <- list(tissue_type = selection)
    expect_false(analysis_hic_enabled(config))
    expect_equal(nrow(load_analysis_hic(config)$PO), 0L)
    expect_equal(nrow(load_analysis_hic(config, all_tissues = TRUE)$PP), 0L)
  }
  config$regulatory$hic <- NULL
  expect_false(analysis_hic_enabled(config))
  config$regulatory$hic <- list(tissue_type = "Gastric")
  local_mocked_bindings(get_hic = function(tissue_type, cache_dir) tissue_type)
  expect_equal(load_analysis_hic(config), "Gastric")
  expect_null(load_analysis_hic(config, all_tissues = TRUE))
  config$regulatory$hic$enabled <- FALSE
  expect_false(analysis_hic_enabled(config))
})

test_that("JEME-only annotations, classification and networks need no HiC resources", {
  skip_if(Sys.which("bedtools") == "")
  config <- write_network_test_config(available = "[network_1, network_2]", selected = "[network_1, network_2]")
  config$aracne <- list(enabled = FALSE)
  config$regulatory$jeme <- list(key_column = "tiss", value = "E094")
  config$regulatory$hic <- list(tissue_type = NULL)
  config$regulatory$promoter_harmonization <- list(enabled = TRUE, gtf_file = "gtf.rds")
  config$reference$ncbi_feature_table <- "ncbi.txt"
  config$networks$split_multi_promoter_edges <- TRUE
  config$networks$deduplicate_jeme_edges_before_harmonization <- TRUE
  config$enhancer_classification <- list(enabled = TRUE, source_analysis = ".", enhancer_flank_bp = 0L,
    promoter_flank_bp = 0L, promoter_window = "strand_aware_tss",
    class1_label = "enh_class1", class2_label = "enh_class2", class3_label = "enh_class3")
  config$ldak$promoter_interval_mode <- "strand_aware_tss"
  config$ldak$annotation_tissues <- list(jeme = "ALL", hic = character())
  readr::write_lines("fixture", file.path(config$project_root, "ncbi.txt"))
  saveRDS(tibble::tibble(gene_id = "ENSG1", gene_name = "A"), file.path(config$project_root, "gtf.rds"))
  jeme <- tibble::tibble(enhancer = "chr1:100-110", CHR = "1", START = 100L, END = 110L,
    promoterFull = "A", promoter = "A", ENSG = "ENSG1", CHR2 = "1", location = 500L,
    desc1 = "Gastric", nfile = 92L)
  local_mocked_bindings(
    get_jeme = function(...) jeme,
    get_hic = function(...) stop("HiC must not be read"),
    NCBI_feature_table = function(...) tibble::tibble(gene_symbol = "A", gene_name = "A",
      CHR = 1L, START = 500L, END = 550L, strand = "+"),
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
  expect_equal(details$source, "JEME")
  for (kind in c("enhancer", "promoter")) {
    readr::write_tsv(tibble::tibble(gene = if (kind == "enhancer") "chr1:100-110" else "A",
      cohort = "test", LRT_P_Perm = 1e-6), file.path(paths[["ldak_summary"]], paste0(kind, "_ldak.tsv.gz")))
  }
  expect_no_error(build_hhotnet_networks(config, dry_run = FALSE))
  edges <- readr::read_tsv(file.path(paths[["hhotnet_data"]], "regulatory_edges.tsv"), show_col_types = FALSE)
  expect_equal(edges$source, "JEME")
  expect_equal(nrow(edges), 1L)
})
