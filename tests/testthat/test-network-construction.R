test_that("network construction dry run respects one-or-more requested networks", {
  config <- write_network_test_config(
    available = "[network_1, network_2, network_3, network_4]",
    selected = "[network_1, network_2, network_3, network_4]"
  )
  single <- build_hhotnet_networks(config, networks = "network_4", dry_run = TRUE)
  multiple <- build_hhotnet_networks(config, networks = c("network_2", "network_4"), dry_run = TRUE)
  expect_named(single, "network_4")
  expect_true(grepl("network_4_edge_list.tsv$", unname(single)))
  expect_named(multiple, c("network_2", "network_4"))
  expect_error(build_hhotnet_networks(config, networks = "network_5", dry_run = TRUE), "not available")
})

test_that("score transformation and graph construction retain validated endpoints", {
  expect_equal(score_from_pvalue(c(0.01, 0)), c(2, 3))
  expect_error(score_from_pvalue(1.1), "must lie")
  nodes <- data.frame(gene = c("A", "B", "C"), score = c(4, 2, 1), node_type = "test")
  edges <- data.frame(Regulator = c("A", "B"), Target = c("B", "C"), source = "JEME", weight = NA_real_)
  graph <- make_network_graph(edges, nodes)
  small <- subset_network_components(graph, threshold = 4)
  expect_equal(igraph::vcount(small), 3)
  expect_equal(igraph::ecount(small), 2)
  expect_error(make_network_graph(data.frame(Regulator = "A", Target = "Z", source = "JEME"), nodes), "missing node")
})

test_that("gene-node extraction supports a disabled empty ARACNe source", {
  jeme_nodes <- list(promoters = data.frame(promoter = "A"))
  hic_nodes <- list(promoters = data.frame(Promoter = "B"))
  empty_aracne <- data.frame(Regulator = character(), Target = character())
  result <- hotnetR2:::extract_gene_nodes(jeme_nodes, hic_nodes, empty_aracne)
  expect_setequal(result$gene, c("A", "B"))
  expect_false(any(result$sources == "ARACNe"))
})

test_that("self-loop policy applies consistently before graph construction", {
  edges <- data.frame(
    Regulator = c("A", "A", "B"),
    Target = c("A", "B", "B"),
    source = c("ARACNe", "JEME", "HiC_PP")
  )
  expect_identical(apply_network_self_loop_policy(edges, TRUE), edges)
  without_loops <- apply_network_self_loop_policy(edges, FALSE)
  expect_equal(nrow(without_loops), 1)
  expect_equal(without_loops$source, "JEME")

  nodes <- data.frame(gene = c("A", "B"))
  kept <- filter_network_regulatory_edges(edges, nodes, "current", TRUE)$kept
  dropped <- filter_network_regulatory_edges(edges, nodes, "current", FALSE)$kept
  expect_equal(nrow(kept), 3)
  expect_equal(nrow(dropped), 1)
})

test_that("ARACNe harmonization maps Ensembl and aliases conservatively", {
  mapping_dir <- tempfile("aracne-mapping-")
  dir.create(mapping_dir)
  saveRDS(
    data.frame(
      gene_id = c("ENSG000001.2", "ENSG000002.1", "ENSG000003.1"),
      gene_name = c("GENE1", "GENE2", "GENE3")
    ),
    file.path(mapping_dir, "gtf.rds")
  )
  saveRDS(
    data.frame(
      ensg = c("ENSG000001", "ENSG000002", "ENSG000009", "ENSG000009"),
      gene = c("GENE1", "GENE2", "GENE9A", "GENE9B"),
      Hsym = c("OLD1", "GENE2", "OLD9", "OLD9")
    ),
    file.path(mapping_dir, "aliases.rds")
  )
  edges <- data.frame(
    Regulator = c("ENSG000001.2", "OLD1", "GENE2", "ENSG000009", "UNKNOWN"),
    Target = c("GENE2", "GENE2", "OLD1", "GENE2", "GENE3"),
    MI = 1:5,
    pvalue = rep(1e-8, 5)
  )
  result <- harmonize_aracne_endpoints(
    edges,
    list(gtf_file = "gtf.rds", alias_file = "aliases.rds"),
    mapping_dir
  )
  expect_equal(nrow(result$data), 4)
  expect_equal(result$data$Regulator[1:2], c("GENE1", "GENE1"))
  expect_equal(result$data$Regulator_mapping_status[1:2],
               c("ensembl_mapped", "alias_mapped"))
  expect_equal(result$data$Regulator[4], "UNKNOWN")
  expect_equal(result$data$Regulator_mapping_status[4],
               "unmatched_symbol_retained")
  expect_false(any(result$data$Regulator_original == "ENSG000009"))
  expect_equal(result$edge_summary$dropped_edges, 1)
  expect_true(all(c("Regulator_original", "Target_original") %in%
                    names(result$data)))
})

test_that("ARACNe harmonization is enabled by default in configuration", {
  config <- write_network_test_config()
  settings <- default_aracne_harmonization(config$aracne$harmonization)
  expect_true(settings$enabled)
  expect_equal(settings$ambiguous, "drop")
  expect_equal(settings$unmapped_ensembl, "drop")
  expect_equal(settings$unmatched_symbols, "keep")
})

test_that("MHC exclusion removes every node type with any inclusive overlap", {
  nodes <- data.frame(
    gene = c(
      "inside_enhancer", "left_overlap", "right_overlap",
      "left_boundary", "right_boundary", "outside_left",
      "outside_right", "other_chr", "duplicate", "duplicate"
    ),
    node_type = c(
      "enhancer", "promoter_gene", "other",
      "enhancer", "promoter_gene", "enhancer",
      "promoter_gene", "other", "enhancer", "enhancer"
    ),
    CHR = c(6, "chr6", 6, 6, 6, 6, 6, 5, 6, 5),
    START = c(
      26000000, 24900000, 33900000, 24000000, 34000000,
      24000000, 34000001, 26000000, 26000000, 100
    ),
    END = c(
      27000000, 25000000, 34100000, 25000000, 35000000,
      24999999, 35000000, 27000000, 27000000, 200
    )
  )
  result <- exclude_network_mhc(
    nodes,
    list(
      enabled = TRUE, genome_build = "hg19", chr = 6,
      start = 25000000, end = 34000000
    )
  )
  expect_setequal(
    result$excluded$gene,
    c(
      "inside_enhancer", "left_overlap", "right_overlap",
      "left_boundary", "right_boundary", "duplicate"
    )
  )
  expect_setequal(
    result$data$gene,
    c("outside_left", "outside_right", "other_chr")
  )
  expect_equal(result$region$start, 25000000)
  expect_equal(result$region$end, 34000000)
})

test_that("MHC exclusion fails closed for invalid coordinates", {
  nodes <- data.frame(
    gene = c("valid", "unknown"),
    node_type = c("enhancer", "promoter_gene"),
    CHR = c(6, 6),
    START = c(25000000, NA),
    END = c(25000001, 26000000)
  )
  expect_error(
    exclude_network_mhc(nodes, list(enabled = TRUE)),
    "cannot guarantee complete removal"
  )
})

test_that("MHC exclusion configuration is fixed to the conservative hg19 region", {
  config <- write_network_test_config()
  config$regulatory$mhc_exclusion <- list(
    enabled = TRUE,
    genome_build = "hg19",
    chr = 6,
    start = 25000000,
    end = 34000000
  )
  expect_silent(validate_analysis_config(config))
  expect_silent(build_hhotnet_networks(config, dry_run = TRUE))

  narrower <- config
  narrower$regulatory$mhc_exclusion$start <- 28477797
  expect_error(validate_analysis_config(narrower), "fixed at 25000000")

  selective <- config
  selective$regulatory$mhc_exclusion$promoters <- TRUE
  expect_error(validate_analysis_config(selective), "applies to all node types")

  wrong_build <- config
  wrong_build$reference$genome_build <- "GRCh38"
  wrong_build$reference$annotation_genome_build <- "GRCh38"
  expect_error(validate_analysis_config(wrong_build), "requires hg19/GRCh37")
})

test_that("MHC-excluded nodes cannot remain in regulatory or ARACNe edges", {
  nodes <- data.frame(
    gene = c("mhc_enhancer", "mhc_promoter", "kept_a", "kept_b"),
    node_type = c("enhancer", "promoter_gene", "enhancer", "promoter_gene"),
    CHR = c(6, 6, 5, 5),
    START = c(25000000, 33900000, 10, 30),
    END = c(25000001, 34000000, 20, 40)
  )
  kept_nodes <- exclude_network_mhc(nodes, list(enabled = TRUE))$data
  edges <- data.frame(
    Regulator = c("mhc_enhancer", "kept_a", "kept_a"),
    Target = c("kept_a", "mhc_promoter", "kept_b"),
    source = c("JEME", "ARACNe", "HiC_PP")
  )
  filtered <- filter_network_regulatory_edges(edges, kept_nodes, "current")
  expect_equal(nrow(filtered$kept), 1)
  expect_equal(filtered$kept$Regulator, "kept_a")
  expect_equal(filtered$kept$Target, "kept_b")
  expect_equal(nrow(filtered$dropped), 2)
})

test_that("class2 and class3 JEME and HiC-PO edges are removed before network construction", {
  config <- write_network_test_config()
  config$enhancer_classification <- list(
    enabled = TRUE, class1_label = "enh_class1", class2_label = "enh_class2",
    class3_label = "enh_class3"
  )
  annotations <- analysis_paths(config)[["ldak_annotations"]]
  dir.create(annotations, recursive = TRUE)
  readr::write_tsv(
    data.frame(
      name = c("enh_keep", "enh_drop_contained", "enh_drop_partial"),
      enhancer_class = c("enh_class1", "enh_class2", "enh_class3")
    ),
    file.path(annotations, "enhancer_classification.tsv.gz")
  )
  data <- list(
    jeme = list(data.frame(
      enhancer = c("enh_keep", "enh_drop_contained", "enh_drop_partial"),
      promoter = c("A", "B", "C"),
      CHR = 1L, START = c(10L, 30L, 50L), END = c(20L, 40L, 60L),
      promoterFull = c("ENSG1$A$chr1$100$+", "ENSG2$B$chr1$200$+", "ENSG3$C$chr1$300$+"),
      ENSG = c("ENSG1", "ENSG2", "ENSG3"), CHR2 = "chr1", location = c(100L, 200L, 300L),
      desc1 = "fixture"
    )),
    hic = list(
      PO = data.frame(
        Interacting_fragment = c("enh_keep", "enh_drop_contained", "enh_drop_partial"),
        Promoter = c("D", "E", "F"),
        CHR = 1L, START = c(10L, 30L, 50L), END = c(20L, 40L, 60L),
        Tissue_type = "fixture"
      ),
      PP = data.frame(Promoter1 = "A", Promoter2 = "B", Tissue_type = "fixture")
    ),
    aracne = data.frame(Regulator = "X", Target = "Y"),
    nodes = list()
  )
  result <- classify_network_enhancer_edges(data, config)
  expect_equal(nrow(result$data$jeme[[1]]), 1L)
  expect_equal(result$data$jeme[[1]]$enhancer, "enh_keep")
  expect_equal(nrow(result$data$hic$PO), 1L)
  expect_equal(nrow(result$data$hic$PP), 1L)
  expect_equal(sum(result$audit$edge_action == "DROP"), 4L)
  expect_setequal(
    unique(result$audit$enhancer_class[result$audit$edge_action == "DROP"]),
    c("enh_class2", "enh_class3")
  )
})
