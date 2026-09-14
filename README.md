# hotnetR2

A standalone development successor to hotnetR combining dataset preparation with
the v2_analysis2 and v3_analysis2 LDAK / hierarchical HotNet workflows.

**Status: implementation in progress (0.0.3); clean-machine reproduction is
not yet certified.** No existing analysis is modified. When implementations
conflict, v2/v3 behavior takes precedence over earlier analysis2 behavior.

## Installed-package workflow

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
