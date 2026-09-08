# hotnetR2

A standalone development successor to hotnetR combining dataset preparation with
the v2_analysis2 and v3_analysis2 LDAK / hierarchical HotNet workflows.

**Status: implementation in progress (0.0.1); clean-machine reproduction is
not yet certified.** No existing analysis is modified. When implementations
conflict, v2/v3 behavior takes precedence over earlier analysis2 behavior.

## Installed-package workflow

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
