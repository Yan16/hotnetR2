# hotnetR2 — current specification and decision ledger

Audit date: 2026-09-28. Runtime baseline: **0.0.7, a722b47** (clean source at
audit start); current source version: **0.2.4**. This document describes
accepted/current behavior, not a proposal to change science. Documentation and
version commits do not change the frozen 0.0.7 runtime baseline.
Read with `KNOWN_ISSUES.md`, `CHANGELOG.md`, and the sibling
`../buas_paper_2026_v2/PROJECT_SPEC.md`. `NEWS.md` remains the version history.

## Authority and boundaries

Latest explicit user decisions define intent; current executable code and YAML
define implementation; checks/results establish verification. Disagreement is
an issue, not permission to silently fix a parameter. v2/v3 decisions supersede
earlier analysis2 where they conflict. v6/v7 are current reproduction targets,
not universal package defaults. Do not overwrite their datasets or results.

Three layers are separate: package/reference preparation; external GWAS/ARACNe
preprocessing; project-specific orchestration. Package functions use explicit
configuration, not the current directory or original author's analysis folders
as implicit scientific inputs. Some portability gaps remain (see issues).

Dataset assemblies, coordinate systems and conversion-boundary decisions are
recorded in `docs/dataset_coordinate_reference.md`. In particular, GWAS and
PLINK single-base variant positions are 1-based identifiers, whereas LDAK
`.loc` intervals follow its 0-start, half-open contract. Coordinate metadata is
field-specific: JEME `location` is a 1-based hg19 TSS, JEME enhancer bounds are
already BED, and active v6/v7 JEME/HiC promoter bounds are obtained by matching
harmonized symbols to the 1-based-inclusive NCBI GRCh37 feature table.

## Current architecture

```text
setup_data / curated reference builder / provision_resources
                      ↓ frozen input cache
read_analysis_config → setup_analysis → annotations → classification
     → LDAK → score summaries → regulatory + interaction networks
     → HHN similarity
        ├─ permutation mode: observed + null hierarchies → selected cut → exports
        └─ manual_delta: stop; exports builds observed hierarchy → supplied cuts
```

Source map (paths relative to package):

| Responsibility | Implementation |
| --- | --- |
| Download/cache readers | `R/download_JEME.R`, `download_HiC.R`, `download_GTEx.R`, `JEME_access.R`, `data.R` |
| Setup and harmonization resources | `R/setup.R`, `harmonization_setup.R`, `inst/scripts/setup_data.R`, `inst/scripts/harmonization/` |
| Configuration, validation, paths | `R/workflow_config.R` |
| User-facing schema reference | `vignettes/analysis-yaml-reference.Rmd` |
| Portable bundles and historical templates | `R/resources.R`, `inst/extdata/v2.yaml`, `v3.yaml` |
| Annotation coordinates and source harmonization | `R/workflow_ldak_annotations.R`, other harmonization helpers in `R/` |
| Source enablement/filtering | `R/workflow_sources.R`, `hic_jeme_overlap.R` |
| Enhancer classes and retention | `R/classification.R`, `enhancer_retention.R` |
| LDAK jobs and summaries | `R/workflow_ldak_jobs.R`, `workflow_ldak_summary.R` |
| Scored graphs | `R/workflow_network_construction.R` |
| STRINGdb cache/setup and augmentation | `R/stringdb_workflow.R` |
| Execution/reuse | `R/execution.R`, `workflow_hhotnet_jobs.R` |
| Cluster parsing/export | `R/workflow_hhotnet_results.R`, `workflow_hhotnet_completion.R` |
| Manual cuts | `R/hhotnet_delta.R`, `inst/scripts/cut_hhotnet_delta.py` |

Preserve public names/arguments when semantics remain the same. Use tidyverse
with namespace-qualified functions; base R/stats are acceptable. Roxygen is
required for public functions. Validate external inputs and configuration;
avoid duplicate defensive checks on internal values already guaranteed by the
producer. Keep original Stata/R/Rmd preprocessing programs as external assets,
not code executed during load/install. No hidden network access in installation.

## Decision ledger: reference inputs and harmonization

Every row is current unless marked historical. Evidence is code inspection and
release tests named in `tests/testthat/`, supplemented by the sibling test repo.

| ID | Accepted behavior / scope | Evidence; superseded behavior |
| --- | --- | --- |
| P01 | `setup_data(cache_dir=NULL, overwrite=FALSE, create_harmonization=FALSE, ...)` is explicit setup; installation does not refresh caches. | `setup.R`, `inst/scripts/setup_data.R`; `setup_hotnet_data` is deprecated. |
| P02 | JEME lasso is the target profile; download API also supports elasticnet. HiC PO/PP tables are separate caches. | download/readers; do not collapse source-specific edge types. |
| P03 | JEME gene mapping is built from target ENSG across all tissues, then applied to annotations and selected-tissue networks. HiC maps symbols/aliases. | annotation/network harmonization; supersedes direct JEME-symbol-only mapping. |
| P04 | September 2026 final curated workbooks define the packaged alias refresh. C11orf48 → LBHD1; MEGT1 → LY6G6D. | external builder and manual-resolution audit; no global C11orf98 or fusion-locus replacement. |
| P05 | New tables carry `mapping_direction=raw_to_approved`; unmarked existing caches retain legacy interpretation. Do not strip attributes or silently replace old caches. | harmonization helpers; critical distinction between fresh defaults and v6/v7's frozen cache. |
| P06 | `data(alias_link)` loads object `alias`; `data(alias_link_nodup)` loads `alias_nodup`. Object name is not necessarily the data-set name. | packaged data/Rd; `data(alias_nodup)` is not the documented loading contract. |
| P07 | HiC with empty/omitted tissue or `enabled=FALSE` is disabled in workflow loaders; public `get_hic(NULL)` still means all tissues. | `workflow_sources.R`, JEME-only tests; prevents formerly mandatory HiC selector errors. |
| P08 | `hic.edge_types` omitted means PO+PP. Accept PO, PP, or both; reject empty/duplicate/unknown entries. Unselected cache is not read. | source loader and PP-only tests. |

### Reference creation versus restoration

`setup_data` invokes JEME, HiC, NCBI, GTEx/GENCODE, HGNC downloaders. Source
locations are embedded in their modules (URLs not revalidated online in this
audit). Freeze retrieved bytes; a live HGNC endpoint is not a version lock.

- JEME: YipLab `encoderoadmap_lasso.zip`; processed per-file RDS under
  `.cache/JEME/lasso/encoderoadmap_lasso.<n>.rds`; tissue-code metadata is also
  packaged. `tiss=E094` is a metadata selection, not an RDS filename of E094.
- HiC: Jung supplementary Excel tables MOESM3/MOESM4, DOI
  `10.1038/s41588-019-0494-8`; `.cache/HiC/Jung_HiC_P-O.rds` and `P-P.rds`.
- NCBI: GRCh37.p13 feature table, release path `GCF_000001405.25-RS_2024_09`;
  `.cache/LDAK/GCF_000001405.25_GRCh37.p13_feature_table.txt.gz`.
- GENCODE v39/GRCh38 (GTEx v10 reference) supplies identifier harmonization,
  **not GRCh37 promoter coordinates**. Cleaned `.cache/gencode.v39.GRCh38.genes_nodup.rds`.
- `create_harmonization_references(cache_dir, overwrite=FALSE, alias_file=NULL)`
  builds GTF and alias caches. `create_alias_references()` can build aliases only.
  Full reference creation requires rtracklayer; explicit legacy alias workbook
  preserves its historical interpretation.

External alias rebuild, from an installed package:

```r
source(system.file("scripts/harmonization/build_alias_data.R", package="hotnetR2"))
build_alias_data(
  source_dir=system.file("extdata/harmonization_2026_0909", package="hotnetR2"),
  legacy_file=system.file("extdata/gtex_hugo_merged_by_ens_SORTens_FINALv1.xlsx", package="hotnetR2"),
  output_dir="rebuilt_aliases")
```

Final workbooks, original programs, source hashes and resolution audit are in
`inst/extdata/harmonization_2026_0909/`. Original Stata scripts need external
inputs and are not certified equivalent to the original R scripts. The portable
builder consumes already-curated final workbooks; further conflicts stop it.

`capture_resources(root, files, file=NULL)` records SHA256; `provision_resources`
accepts a local bundle or URL manifest, verifies hashes and refuses different
existing content unless overwrite is explicit. Relative resource paths cannot
contain parent traversal. These are restoration tools, not provenance evidence
for how an undocumented dataset was originally made.

## Decision ledger: intervals and enhancer policy

| ID | Current rule | Evidence / critical qualification |
| --- | --- | --- |
| P09 | All `.loc` and bedtools-facing intervals use 0-based, half-open BED coordinates. Convert each 1-based source exactly once and never decrement an already-BED enhancer start. | Coordinate regression fixtures cover both strands, gene bodies, zero starts and touching boundaries. Supersedes frozen legacy arithmetic in 0.1.1. |
| P10 | NCBI TSS conversion uses `TSS=START` for `+` and `TSS=END` for `-`, then writes `[TSS-1,TSS)`. NCBI gene bodies write `[START-1,END)`. Original biological bounds and strand remain in detail tables. | `workflow_ldak_annotations.R`; invalid strand fails. |
| P11 | Classifier expands stored BED intervals directly: start `max(0,START-flank)`, end `END+flank`. Same chromosome only; touching half-open boundaries do not overlap. | `classification.R`; no second start conversion. |
| P12 | enh_class1: no overlap with any promoter window; class2: entire expanded enhancer contained in at least one window; class3: other overlap. | class tests; supersedes earlier two-class categorization. Not restricted to JEME target pairs. |
| P13 | Default drop sources: TSS → HiC_PO only; gene body → JEME and HiC_PO. Explicit list overrides; empty retains all. A mixed-source annotation is retained if any source retains it, but edge filtering remains source-specific. | `enhancer_retention.R`; replaces universal class1-only TSS behavior. |
| P14 | Optional `hic.jeme_overlap.enabled=TRUE`: compare HiC PO BED intervals to the **selected** JEME BED tissue union, default flank 1000 on each set; any positive-width overlap removes all PO contacts for that enhancer. Boundary-only contact is not overlap. PP and JEME are unchanged. | `hic_jeme_overlap.R`; works before annotations AND networks, even if annotation universe is all tissues. |
| P15 | HiC–JEME class1/2 is separate from enhancer–promoter class1/2/3. Passing first filter does not bypass TSS or score filters. | source/classification integration; v6/v7 enable both. |

`.loc` is headerless `name, CHR, START, END`, unexpanded. LDAK receives
`--gene-buffer` separately. Do not feed already-expanded overlap-audit intervals
to LDAK. Chromosomes normalize X/Y to 23/24; NCBI reader excludes mitochondrial
genes, sorts and keeps one row per symbol, not a chosen longest transcript.
Annotation validation checks assembly declarations, names, numeric intervals,
duplicates and reference bounds. Changing duplicate selection changes coordinates.

## Decision ledger: scores, networks, and HHN

| ID | Rule | Evidence |
| --- | --- | --- |
| P16 | Node score = `-log10(p)` from configured column. With null cap, positive infinity becomes maximum finite score +1 (or 1 if none). Invalid p outside [0,1] fails. | `score_from_pvalue`; root `top_enhancer.R` is a report and deliberately leaves Inf, not the HHN scoring implementation. |
| P17 | Network1/2 retain entire undirected regulatory connected components with maximum score >= threshold. Weak nodes can remain. | `subset_network_components`; thresholds are not per-node pruning or p-value filters on every edge. |
| P18 | Network3/4 add ARACNe edges with MI >= threshold between nodes already in Network1/2, respectively; they do not introduce new ARACNe-only endpoints. | `add_network_aracne_edges`; initial ARACNe reader applies MI>=0 and pvalue<=1e-6 in target profiles. |
| P19 | Self-loops default retained; legacy-input ordering and source-specific dedup/splitting in v6/v7 are part of reproduction. `preserve_duplicate_edges_during_aracne` does not mean no upstream deduplication. | `workflow_network_construction.R`; node IDs also use punctuation-normalized matching. |
| P20 | Permutation mode is default. Local observed index is 0; score permutations use base seed+i. Current local profile base 0 gives seeds 1..100. | generated HHN script; Apptainer seed equivalence not established. |
| P21 | manual_delta is opt-in, local_python only. `hotnet` stops after similarity; `exports` creates/reuses observed score-weighted hierarchy, calls HHN `cut_hierarchy`, no nulls. | `hhotnet_delta.R`, adapter, synthetic tests. |
| P22 | Delta default vector c(.05,.1,.2,.5), positive finite; duplicate values collapse; filename collisions at 15-digit representation fail. Cut includes heights >= delta. | native HHN function; larger delta generally yields smaller clusters, not guaranteed biological significance. |
| P23 | Raw manual TSV includes all clusters/singletons; CX2 min_cluster_size defaults 2, only within-cluster original edges. Empty graph is valid. Standard permutation export instead includes connected components containing multi-node clusters. | result exporter versus delta exporter; do NOT unify these selection rules accidentally. |
| P24 | HHN hierarchy construction restricts to a largest SCC when needed. Not every input-network node must be in reported clusters. | local HHN `construct_hierarchy.py`; v6 network4: 2391 graph nodes, 1732 hierarchy members. |
| P25 | Interaction backends are explicit. Omitted `interaction_network.type` preserves legacy ARACNe selection. With `type: stringdb`, Network3/4 augment Network1/2 using cached human STRING v12 interactions with `combined_score >= 700` for the current profile. Only already-present `promoter_gene` nodes may map; no enhancers or new endpoints are introduced. Multiple genes for one STRING ID resolve by shortest name then lexical order and are audited. STRING loops, duplicate gene pairs and pairs already present in the base graph are not added, so existing regulatory edges and sources remain intact. | `stringdb_workflow.R`, `test-stringdb-workflow.R`; cache download is an explicit setup action, never an analysis side effect. |
| P26 | Documentation and new analysis templates use the corrected E092/Gastric TSS profile: NCBI GRCh37 strand-aware TSS anchors, 0-based half-open intervals, no enhancer-class exclusion, all-tissue JEME ENSG mapping, selected-tissue JEME/HiC edges and current curated aliases. v6/v7 remain references for four-network orchestration and permutation/manual-delta modes, not exact corrected-coordinate or refreshed-alias baselines. | `vignettes/`, `inst/extdata/vignette_examples/`, example tests. Gene-body/gene-window plus enhancer-class exclusion is a separately named historical alternative; original analysis2 and early intermediates are obsolete templates. |
| P27 | LDAK summarization also creates `enhancer_ldak_long.tsv.gz` by left-expanding the standardized enhancer table over distinct all-tissue JEME enhancer--target associations. It preserves raw JEME `promoterFull`, `ENSG`, and `original_promoter`, adds the configured-harmonized `promoter` plus tissue metadata, and retains non-JEME enhancers with missing target fields. The standard enhancer table stores sorted distinct raw promoter labels collapsed with `;`. | `create_enhancer_ldak_long_summary()`, summary tests and real-cache test. The updater also upgrades an existing standard enhancer summary without rerunning or changing LDAK statistics. |
| P28 | `Gene_Name` is the sole public node identifier in standardized enhancer and promoter LDAK summaries. The redundant duplicate `gene` column is not written. Current network construction reads `Gene_Name` and only renames it to the separate internal network-node `gene` schema after summary ingestion. | `summarize_one_ldak_result()`, `create_enhancer_ldak_long_summary()`, `filter_analysis_network_ldak()` and regression tests. A legacy enhancer summary is accepted only when `gene` exactly agrees with `Gene_Name`; the updater then removes `gene`. |

Key API contract:

```r
config <- hotnetR2::read_analysis_config("analysis.yaml")
hotnetR2::run_analysis(config) # annotations, classification, ldak, scores, networks, hotnet, exports
hotnetR2::run_analysis(config, stages=c("networks", "hotnet", "exports"))
hotnetR2::extract_hhotnet_delta_clusters(config,
  deltas=c(.05,.1,.2,.5), networks=NULL, min_cluster_size=2L)
```

`networks=NULL` uses YAML selection. Selected network identifiers are validated;
Network3/4 internally require Network1/2 bases. `generate_hhotnet_job(...,
dry_run=TRUE)` returns path without writing; explicit FALSE writes an executable
runner. `prepare_ldak_annotations`/`build_hhotnet_networks` also have dry-run APIs;
consult installed Rd before direct execution. `initialize_analysis(profile)`
currently accepts only historical `v2` or `v3`, NOT v6/v7. Use audited v6/v7 YAML
snapshots for new target profiles rather than assuming initializer equivalence.

## Defaults are not profile settings

| Setting | Package fallback / acceptance | v6/v7 |
| --- | --- | --- |
| promoter_interval_mode | gene_body; strand_aware_tss supported | strand_aware_tss |
| HiC tissue / edge types | empty disables; absent types means PO+PP | Gastric; PO+PP |
| HiC–JEME exclusion | opt-in; flanks default 1000 | enabled, 1000/1000 |
| Interaction backend | omitted means legacy ARACNe selection | ARACNe for v6/v7; STRINGdb is an explicit new profile |
| STRINGdb release / species / score | 12.0 / 9606 / 700 when STRINGdb is selected | not used by v6/v7 |
| MHC exclusion | FALSE if omitted | FALSE |
| HHN analysis mode | permutation | v6 permutation, v7 manual_delta |
| HHN execution mode | apptainer | local_python |
| HHN seed / python | 0 / python3 | 0 / python3 |
| Annotation tissue universe | implemented target mode ALL | ALL, selected tissues for network/filter reference |
| Canonical documentation profile | corrected strand-aware TSS; enhancer classification disabled | frozen v6/v7 are orchestration references only |

Other required numerical settings are explicit in target YAML. Do not invent
universal defaults for values merely present in an example template.

## I/O, validation, and execution invariants

- Analysis-local `LDAK/{annotations,scripts,results,summary}` and
  `hHotnet/{data,scripts,intermediate,results,summary,cytoscape}`; metadata/logs
  local. Shared datasets are input-only in migration. Absolute author paths are
  not portable configuration.
- Standard summaries: `enhancer_ldak.tsv.gz`, `promoter_ldak.tsv.gz`; fields
  include raw LDAK statistics, canonical node identifier `Gene_Name`, node_type,
  cohort, flank and annotations. The former duplicate `gene` field is absent.
  `enhancer_ldak.tsv.gz` additionally contains `original_promoter`, with sorted
  distinct raw JEME promoter labels collapsed by `;` per enhancer.
  `enhancer_ldak_long.tsv.gz` repeats the enhancer fields for each unique JEME
  `promoterFull`/`ENSG`/`original_promoter`/harmonized-promoter/tissue
  association and retains unmatched enhancers with missing JEME fields. Network
  data: headerless contiguous 1-based index and two-column
  integer edge list; named score TSV; full node/edge tables retain attributes.
- STRINGdb profiles additionally write `stringdb_mapping.tsv`,
  `stringdb_mapping_conflicts.tsv`, `stringdb_summary.tsv` and
  `stringdb_edges.tsv`. The three pinned STRING cache files are validated as
  inputs and included in analysis signatures outside the LDAK-only signature.
- Cluster TSV: comments start `#`; each nonblank noncomment line is a cluster,
  tab-separated node IDs. Reject duplicate membership. Cluster numbers are
  output-local, not stable biological identifiers across runs/deltas.
- Standard export prefix is `config$analysis$name_network_N_hotnet_components`.
  Folder-derived names are a **project launcher** rule, not package behavior.
- Manual outputs reside in separate `manual_delta` subdirectories; filenames
  contain `_delta0.1`, etc.; latest-call manifest overwrites previous manifest,
  while older delta files remain. Do not interpret the directory as one grid.
- Content fingerprints gate reuse; old/unstamped changed result directories
  are archived in `superseded-*`. LDAK refuses linked result write targets.
  An interrupted run is not certified partial resume. External exit 0 alone
  does not prove LDAK success; expected artifacts are checked.
- `check_analysis_inputs(write_report=FALSE)` is read-only presence validation,
  not a complete cache/environment/provenance verification.
- `validate_hhotnet_completion` checks matrix existence in manual mode; it does
  not certify manual grid completion. Verify manifest and files separately.

## Reproduction and development acceptance

Freeze inputs, source, tools and output hashes before refactoring; use the
sibling project specification's comparison scope. Exact identity of selected
v2/v3 artifacts in historical reports is NOT evidence of whole-folder or
fresh-machine identity. v2 fresh promoter permutation variability remains open.

Saved package release checks: `../test_hotnetR2/reports/` and `logs/`; 0.0.7
passed metadata-free R CMD check with vignettes and installed-package synthetic
HHN tests. No new package build/check was run for this documentation-only audit.
Future code releases: focused/full tests → Roxygen → metadata-free build/check
→ temporary installed tests → explicit changed-file review → scoped Git commits.
Do not auto-install into the user's library. No release tag unless requested.

Start a new session with this specification, its issues, sibling project spec
and evidence snapshots. First establish characterization tests; then consolidate
configuration/coordinate/filter decisions without changing these invariants.
