# Public API and coding contract

## Preserve names and arguments where semantics match

| Existing API/source | hotnetR2 destination | Decision |
|---|---|---|
| get_jeme(), get_hic() | Same names and existing leading arguments | Add optional explicit cache/resource context at the end. |
| extract_jeme_nodes(), extract_hic_nodes() | Same | Preserve schemas and ordering. |
| NCBI_feature_table(file, keep_classes) | Same | Preserve legacy parsing; add explicit separate interval conversion. |
| harmonize_promoter_by_ensg(), harmonize_hic_gtf() | Compatibility wrappers | Route through one mapping engine with explicit policies. |
| setup_data(), setup_hotnet_data() | Same wrappers | One resource provisioning implementation. |
| download_NBCI_feature_table() | Retained compatibility spelling | Add correctly spelled alias, one implementation. |
| summarize_ldak_results(config_file, ...) in hotnetR | Same legacy API | Preserve recorded full formals. |
| summarize_ldak_results(config, dry_run) in helper | summarize_analysis_ldak_results() | Explicit rename; no argument-shape dispatch. |
| read_analysis_config(), analysis_paths() | Same existing arguments | Add explicit output/cache overrides at the end as needed. |
| prepare_ldak_annotations(), generate_ldak_jobs() | Same | Shared interval and resource contracts. |
| build_hhotnet_networks(), resolve_networks() | Same | Profiles retain selected network behavior. |
| summarize_hhotnet_results(), export_hhotnet_graphs() | Same | Preserve export and cluster contracts. |
| exclude_MHC() duplicate definitions | One function | Record existing signatures/body differences before choosing. |
| filter_network_ldak() in hotnetR versus helper | Preserve public hotnetR API; rename helper internal function | The public table interface and internal data/config interface must not shadow each other. |
| classifier shell script | classify_enhancers(config, ...) | New installed API, bedtools backend. |
| command-line step sequence | run_analysis(config, stages, ...) | New orchestration, invalidation and progress reporting. |
| clean project setup | initialize_analysis(path, profile, ...) | New API for v2/v3 portable templates. |
| input/output auditing | capture_artifacts(), compare_artifacts() | New explicit reproducibility APIs. |

This table defines intent; exact full signatures must be machine-inventoried
before migration. APIs not listed require an explicit migration disposition,
not accidental export loss from copying only selected files.

## Implementation rules

- Use tidyverse for data manipulation and qualified function calls. Avoid
  library() in package functions and broad imports that mask symbols.
- Pure functions transform documented objects; orchestrators perform I/O.
- Validate external YAML, downloaded/user files, tool availability and output
  artifacts. Do not duplicate those checks on every internal call.
- Ambiguity/fallback policies are explicit scientific choices; never infer a
  fallback because a column or file is unexpectedly absent.
- Use explicit join keys and preserve row order where the reference depends
  on it. Compare order as well as membership.
- Keep one owner for ENSG normalization, symbol aliases, interval/flank
  conversion, class precedence, score transformation, and path resolution.
- Export only independently useful workflow operations. Low-level functions
  remain internal and are documented by contracts and focused comments.
- Preserve base/stats use when clear. Qualify other namespaces, including
  utils::, tools::, igraph::, yaml:: and process execution libraries.
- Roxygen for each public function includes title, purpose, every argument,
  return schema, coordinate conventions where relevant, side effects,
  reproducible small examples, and related APIs.
- Tests should exercise meaningful boundaries, invariants and historical
  regressions, not repeat implementation lines as expectations.

## Suggested modules

resources; config; import_jeme; import_hic; import_ncbi; harmonization;
intervals; enhancer_classification; ldak_execution; ldak_summary; scores;
network_construction; hhotnet_execution; cluster_summary; graph_export;
provenance; artifact_comparison; compatibility.

The dependency graph should point from orchestration to pure domain functions
and external adapters, never from core data operations to CLI scripts.
