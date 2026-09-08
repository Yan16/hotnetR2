write_network_test_config <- function(available = "[network_1, network_2, network_3]",
                                      selected = "[network_1, network_3]") {
  config_dir <- tempfile("hotnetR2-config-")
  dir.create(config_dir)
  config_file <- file.path(config_dir, "analysis.yaml")
  writeLines(c(
    "schema_version: 1", "analysis: {name: test}", "project: {root: '.'}",
    "gwas: {dataset: test, summary_file: a, pvalue_file: b, extract_file: c, sample_size: 1}",
    "reference: {ldak_exe: ldak, bfile_prefix: ref, genome_build: GRCh37, annotation_genome_build: hg19}",
    "regulatory: {jeme: {}, hic: {}}", "aracne: {file: aracne}",
    "ldak: {flank_bp: 0, permutations: 0, threads: 1, ignore_weights: true, allow_ambiguous: true}",
    "scores: {pvalue_column: LRT_P_Perm, transform: neg_log10, infinite_score_cap: null}",
    paste0("networks: {available: ", available, ", selected: ", selected,
           ", regulatory_edge_filter: current, node_order: score_descending, ",
           "preserve_self_loops: true, ",
           "preserve_duplicate_edges_during_aracne: false, deduplicate_jeme_edges_before_harmonization: false, ",
           "split_multi_promoter_edges: false, thresholds: {small_regulatory_score: 1, broad_regulatory_score: 1, small_plus_aracne_mi: 0, broad_plus_aracne_mi: 0}}"),
    "hhotnet: {num_permutations: 0, lower_size_bound: 1, cpus: 1}"
  ), config_file)
  read_analysis_config(config_file)
}
