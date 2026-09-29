# hotnetR2

A standalone successor to hotnetR for tissue-informed LDAK and Hierarchical
HotNet analyses. Current examples use corrected E092/Gastric strand-aware TSS
coordinates, no enhancer-class exclusion, current harmonization and explicit
ARACNe or STRINGdb interaction backends.

## Documentation

Start with `vignette("hotnetR2-overview", package = "hotnetR2")`. Companion
guides cover input preprocessing, resource setup, the current corrected
TSS/no-enhancer-class workflow, ARACNe/STRING analysis profiles, HHN execution
and reproduction boundaries. New analyses should copy the validated examples in
`system.file("extdata", "vignette_examples", package = "hotnetR2")`, not a
historical result folder.

**Status: current package release 0.2.1; complete clean-machine reproduction is
not yet certified.** No existing analysis is modified. Historical v2/v3 and
v6/v7 behavior is retained only where the specification explicitly requires it.

## Installed-package workflow

### HiC promoter-promoter-only selection (0.0.5)

```yaml
regulatory:
  jeme:
    key_column: tiss
    value: E094
  hic:
    tissue_type: Gastric
    edge_types: [PP]
```

This retains Gastric JEME enhancer-promoter and HiC promoter-promoter edges.
HiC PO is excluded before harmonization, annotation, classification and network
construction; its cache is not required. A shared JEME enhancer remains a JEME
node. HiC-only PO enhancers and PO-only targets do not enter the annotations.
The existing all-tissue annotation policy still applies to JEME and the enabled
HiC PP source; the selected tissues restrict network edges. Promoter intervals
and enhancer-class retention remain controlled by the usual TSS settings.
An omitted edge_types retains both sources; [PO] and [PO, PP] are also valid.
Empty tissue selection still disables HiC. Start a fresh analysis and calculate
its scores/networks when changing source selection; do not reuse mixed-source
annotations or scores as if they were PP-only results.

### Curated alias update (0.0.4)

The September 2026 JEME/HiC workbooks now update `alias` and `alias_nodup`
(loaded with `data(alias_link)` and `data(alias_link_nodup)`). Explicit curation
maps C11orf48 to LBHD1 and MEGT1 to LY6G6D. Legacy keys absent from the new
workbooks remain available. Original scripts and spreadsheets are preserved
under `system.file("extdata", "harmonization_2026_0909", package="hotnetR2")`;
the portable external builder and its instructions are under
`system.file("scripts", "harmonization", package="hotnetR2")`.

Generate a separate cache without touching existing analyses:

```r
hotnetR2::create_alias_references(cache_dir = "new_alias_cache")
```

Installing the package alone does not replace `.cache/alias_link*.rds`.
Use the new files explicitly in analysis settings when ready to regenerate
affected annotations, LDAK scores and networks. The new references mark their
raw-to-approved mapping direction; older unmarked RDS files keep their original
HGNC-to-GENCODE interpretation. JEME's ENSG-first workflow is unchanged by this
alias-data update; it does not automatically switch to symbol-only matching.

### TSS enhancer retention (0.0.3)

TSS profiles retain class1, class2 and class3 JEME enhancers. Classifications
still describe the same flanked overlaps; only retention changes. HiC-PO remains
class1-only. Gene-body profiles retain their class1-only defaults.

Set `enhancer_classification.drop_non_class1_edge_sources: [HiC_PO]` for this
policy. An explicit `[JEME, HiC_PO]` preserves the previous behavior; `[]`
retains all classes from both sources. If omitted, the default follows the
promoter-window profile. Class-specific `.loc`/detail files remain audits;
`enhancers.loc` and `enhancers.details.tsv.gz` contain the union retained by
either source. Shared JEME/HiC coordinates are scored once, but network edges
are filtered independently by source. Retention does not bypass score thresholds.

When migrating existing runs, regenerate classification and enhancer LDAK scores
before networks/HotNet; class1-only scores cannot represent the newly retained
classes. Use a new enhancer result directory (the v3 template now uses
`enhancer_retained_tss`) and retain old results for provenance.

### JEME-only analyses (0.0.2)

Leave the HiC tissue selection empty to disable HiC:

```yaml
regulatory:
  jeme:
    key_column: tiss
    value: E094
    method: lasso
  hic:
    tissue_type: []
```

`null`, an empty string, or an omitted `hic` section also disable HiC;
`hic.enabled: false` explicitly disables it even with a tissue configured.
No HiC caches or HiC alias mapping are required in this mode. Disabled HiC
contributes no enhancers, promoters, overlap windows, PO edges or PP edges.
All-tissue annotation collection continues for JEME according to the existing
annotation policy; the selected JEME tissue controls network edges.
Regenerate annotations, classification, LDAK scores and networks in a **fresh
output directory** when switching an existing run to JEME-only: previous mixed
annotations and overlap classes are not JEME-only results.

Nonempty HiC settings preserve the previous mixed-source behavior (including
all-HiC annotation collection). The data-reader call `get_hic(NULL)` still
means all tissues; the disabling rule applies to analysis YAML settings.
Custom scripts pinned to package 0.0.1 must be updated to use 0.0.2.

### Run the workflow

```r
config_file <- hotnetR2::initialize_analysis("my-analysis", profile = "v3")
# Supply the frozen inputs and tools at the paths in analysis.yaml.
# Verify a transferred resource bundle before running:
# hotnetR2::provision_resources(manifest, "my-analysis", source_root = bundle)
hotnetR2::run_analysis(config_file)
```

The v2 profile uses extended gene bodies; v3 uses strand-aware TSS windows.
Both retain legacy coordinate arithmetic, classify enhancers with flanks,
keep only class1 regulatory edges, and run Networks 1–2 without ARACNe.
See `vignette("reproduction", package = "hotnetR2")` after installation.

Reusable full-data testing scripts and the actual run commands are retained in
the sibling [`../test_hotnetR2/`](../test_hotnetR2/README.md), including
[`COMMANDS.md`](../test_hotnetR2/COMMANDS.md). Unit tests stay in this package.

Full-data verification has reproduced five selected annotation/harmonization
artifacts and all 12 current network-input TSVs for each profile. Gzip tables
were compared after decompression; uncompressed network TSVs were byte-matched.
These checks reuse frozen raw LDAK outputs. Fresh HotNet runs also reproduced
the checked cluster/export files, allowing only declared path and producer-label
differences. Fresh v3 LDAK results match; v2 promoter permutation-adjusted
p-values differ and remain under investigation. Two historical source-edge audit files are explicitly outside
the current network comparison. See [verification evidence](docs/verification.md).

## Optional HiC–JEME enhancer overlap filter (0.0.6)

```yaml
regulatory:
  jeme:
    key_column: tiss
    value: E094
  hic:
    tissue_type: Gastric
    edge_types: [PO, PP]
    jeme_overlap:
      enabled: true
      hic_flank_bp: 1000
      jeme_flank_bp: 1000
```

The reference is the union of enhancers in the **selected JEME tissues**, not
all JEME tissues. Expand each enhancer at both ends, clip starts to 1, and
intersect on the same chromosome (1-based inclusive source coordinates).
`hic_jeme_class1` has no overlap and retains its PO edges; `hic_jeme_class2`
has any overlap and drops all its PO edges. PP contacts and JEME contacts are
unchanged. The existing enhancer-versus-promoter class filter still applies
afterward; the new class names describe a separate comparison.

Filtering precedes both LDAK annotation construction and network construction.
All-tissue annotation mode still uses only selected JEME tissues as the overlap
reference. Without this setting, behavior is unchanged. bedtools is required.
Under LDAK annotation outputs, `hic_jeme_{annotations,networks}_*.tsv.gz`
records classification (including flank lengths, selection and counts), overlap
pairs, and the expanded selected JEME reference. Coordinates in these audits
are 1-based inclusive. Retained means passing this filter, not guaranteed
survival of subsequent TSS, GWAS or network-score filters.

## Manual HHN delta cuts (0.0.7)

Set `hhotnet.analysis_mode: manual_delta`, `execution_mode: local_python`,
and optionally `delta_values: [0.05, 0.1, 0.2, 0.5]`. The `hotnet` stage stops
after similarity matrices; `exports` builds the observed score-weighted
hierarchy and calls `extract_hhotnet_delta_clusters()` for the grid. No HHN
permutations are run (LDAK permutations are independent). Direct calls accept
a numeric delta vector and optional network vector. Higher delta generally
gives smaller clusters. These are exploratory cuts, not significance tests.

Raw cluster TSVs retain singletons; CX2 and associated tables default to clusters
of at least two nodes and retain only within-cluster edges, not surrounding
connected components. Empty retained sets produce empty CX2s. Files include
`_delta0.05` etc. and live in dedicated `manual_delta` directories. A manifest
under `hHotnet/summary` records the latest grid's paths, counts and method.
Matrix/hierarchy construction may remain costly; repeat cuts reuse the observed
hierarchy via input fingerprints. The HHN source directory must provide
`process_hierarchies.cut_hierarchy` and `construct_hierarchy.py`.

Read in order:

1. [Workflow schematic](docs/workflow.md)
2. [Implementation plan](PLAN.md)
3. [Detailed task list](TASKS.md)
4. [API and coding contract](docs/api.md)
5. [Reproduction contract](docs/reproducibility.md)
6. [Architecture decisions and source findings](docs/decisions.md)
7. [Captured source inventory](docs/inventory/README.md)

The references are `../buas_paper_2026/v2_analysis2` and
`../buas_paper_2026/v3_analysis2`. Original hotnetR source is at `../hotnetR`.
The three repositories are siblings under the same parent directory.
These are development
references, never runtime dependencies of the finished package.

The first release must preserve names and arguments where their meaning is
unchanged, use explicit namespace-qualified tidyverse calls, provide Roxygen
documentation, and reproduce reference artifacts in isolated output trees.
It must install without hotnetR, analysis-local helper source trees, or
access to the original author's directories.

Public reference downloads, user-supplied GWAS inputs, software environments,
and output artifacts have separate versioned manifests. A new machine needs
the frozen inputs or verified source downloads, as well as the R package.

This repository will track reviewable implementation increments. Completed
work and pending acceptance gates are recorded in TASKS.md.
