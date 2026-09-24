# Architecture decisions and inspected source evidence

## ADR 001 — Independent package

Create hotnetR2 as an independently installable package. Original hotnetR and
analysis directories remain baseline sources. Migrate the required dependency
closure rather than depending on hotnetR at runtime or sourcing helper files.
The analysis names requested as v2_analysis/ and v3_analysis/ refer here to the
existing v2_analysis2/ and v3_analysis2/ directories.

## ADR 002 — Freeze behavior before cleanup

Compatibility migration preserves actual generated data, even where a legacy
decision merits correction. Corrections are explicit policies or separately
versioned profiles with their own validation, not invisible refactors.

## ADR 003 — Canonical BED interval policy (implemented in 0.1.1)

The source-coordinate audit is now recorded in
`docs/dataset_coordinate_reference.md`. BEEA `POS`, MVP `pos_hg19`, and PLINK
BIM base-pair positions are confirmed 1-based variant positions. NCBI feature
coordinates are treated as 1-based inclusive at the input boundary. JEME uses
mixed field-level conventions: its Ensembl-derived `location` is a 1-based TSS,
but its enhancer bounds are already BED intervals. HiC promoters have no source
coordinate in the interaction tables; the active workflow assigns their bounds
by gene-symbol/alias matching to NCBI. Jung HiC other-end coordinates describe
the BED-style fragment intervals built from midpoints of adjacent HindIII sites.

Inspected sources:
- hotnetR/R/download_JEME.R: NCBI_feature_table() renames start/end and converts
  to integer without an explicit one-based-to-BED conversion.
- v3 helper R/ldak_annotations.R: convert_promoter_details_to_tss() assumes BED,
  computes plus TSS=START+1, minus TSS=END, writes START=TSS-1 and END=TSS.
- Both classifier shell scripts: expanded start=START-flank, then subtract
  one when writing BED for bedtools, documenting one-based inclusive input.
- Generated LDAK jobs pass .loc directly to LDAK, whose help/output says BED.

Version 0.1.1 resolves the inconsistency: every `.loc` and bedtools-facing
interval is 0-based half-open; NCBI and JEME TSS values and NCBI gene bodies are
converted once; JEME and Jung HiC enhancer/fragment intervals are not converted
again. Touching boundaries are adjacent. Regression fixtures cover both TSS
strands, gene bodies, zero starts and boundary contact. Historical v6/v7 output
directories remain evidence of the superseded arithmetic and must not be
overwritten when corrected analyses are generated.

## ADR 004 — Explicit API collision resolution

hotnetR and helper both export summarize_ldak_results() with different
arguments. Preserve the hotnetR interface and use
summarize_analysis_ldak_results() for config-driven analysis. Duplicate
exclude_MHC() definitions in hotnetR/R/data.R and R/r_utils.R need one owner.

The automated inventory also found filter_network_ldak(): the public hotnetR
function takes (net_regulators, nodes_all, verbose), while the helper's internal
function takes (data, config). Rename the internal orchestration operation;
preserve the public signature and keep the shared filtering decision in one
implementation. The two exclude_MHC() bodies have different hashes despite
matching formals, so source load order currently matters.

## ADR 005 — Explicit resources and immutable identity

Current get_cache_dir() uses here()/working-directory inference. Several
loaders use hotnetR:::, and dataset lookup is hardcoded to package='hotnetR'.
Replace these with package-owned resources and explicit project/cache context.
Historical mappings must not silently change when remote references update.

## ADR 006 — Resume requires fingerprints

The current HotNet runner skips stages based on file presence. New resume
requires matching inputs/config/tool/code fingerprints. LDAK exit status must
be supplemented with output validation. External-input errors are boundary
failures; internal pure transforms should not contain redundant defensive code.

## Source provenance at planning start

- Original hotnetR version: 0.4.6.2.
- Original hotnetR Git HEAD: afde6f826c370bff53e09480ec94cfaf994586cd.
- Containing project HEAD: c6cadde85910b15fa337743ce9685ffe784105f6.
- Hierarchical HotNet checkout HEAD: 448e0d20bb7c8414ab10a4f40cce3cb32944c3a5.
- hotnetR has an existing staged README modification and untracked documents.
  Preserve them; Git HEAD alone does not identify the full working snapshot.
- Record source-tree checksums and installed-package provenance in Phase 0;
  installed-package provenance remains pending. The inventory below records
  the inspected source files and function signatures.

## Repeatable source inventory

tools/capture_source_inventory.R parsed 50 source/metadata files and 323
top-level function definitions without loading either package. Reports are in
docs/inventory/: source_files.tsv, function_signatures.tsv,
duplicate_definitions.tsv and signature_conflicts.tsv. Source file hashes
are raw SHA-256; function-body hashes use digest's R serialization and are
diagnostic within the recorded R environment, not cross-R-version identities.

## Git policy

This new project has its own Git history and lives beside hotnetR and
buas_paper_2026. Its GitHub remote will be configured separately by the owner.
Commit planning first, then focused
implementation increments with task IDs, tests, and evidence. Stage explicit
paths. Do not amend original hotnetR history or include existing workspace
changes. Release tags follow completed package and reproduction gates only.
