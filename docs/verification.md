# Development verification

Reference policy: v2_analysis2 and v3_analysis2 override earlier analysis2.
Original analyses are read-only. Candidate outputs are generated separately.

Run `../test_hotnetR2/scripts/verify_reference_profiles.R PROJECT_ROOT OUTPUT_ROOT` with an installed
hotnetR2 package. The script regenerates annotations and networks; it links raw
reference LDAK result directories as fixed inputs. Never run `run_ldak()` against
those links (the runner rejects them).

## Full-data comparisons obtained

| Profile | Class1 | Class2 | Class3 | Selected annotation files | Network TSVs |
|---|---:|---:|---:|---:|---:|
| v2 | 225571 | 196458 | 24594 | 5/5 identical | 12/12 identical |
| v3 | 424653 | 1951 | 20019 | 5/5 identical | 12/12 identical |

Annotation comparisons cover promoters.loc, enhancers.loc,
enhancers.details.tsv.gz, enhancer_classification.tsv.gz, and
jeme_ensg_harmonization.tsv.gz. Gzip content is compared after decompression;
network TSVs are compared byte-for-byte. The empty ARACNe table is preserved
for schema compatibility, not an enabled ARACNe network.

Historical regulatory_edges_class1_pre_harmonization.tsv and
regulatory_edges_class1.tsv are source-analysis audits, not current network
inputs. Their regeneration remains a separate pending task. No claim of
whole-folder identity follows from these selected comparisons.

Local evidence is under `/Users/li26191/tmp/hotnetR2-reference-runs`.
This development location is not a runtime requirement.

## Pending acceptance gates

- Independently execute LDAK, then compare complete outputs.
- Freeze distributable input provenance and R/Python dependency locks.
- Restore and execute on a clean machine without source-analysis directories.
- Test whole-workflow invalidation; unit tests cover individual fingerprints.

The execution runner reuses only stamped completed outputs. Unstamped or
changed result directories are moved to recoverable `superseded-*` directories;
interrupted stages restart rather than receiving a certified partial resume.

## Fresh HotNet and score verification

The two standardized LDAK summary tables per profile are content-identical
after decompression, using frozen raw LDAK results.

Both profiles were independently executed with 100 HotNet permutations per
network. All four networks produced 101 hierarchy edge/index pairs. Export
validation passed for unique node names, edge endpoints, GraphML readability,
and Cytoscape.js/CX2 counts and attributes.

`../test_hotnetR2/scripts/verify_hotnet_artifacts.R` compares the union of reference and candidate
TSV/GraphML/JSON/CX2 inventories under results, cytoscape, and summary/clusters.
It found 25 files for v2 and 59 for v3: respectively 20 and 54 byte-identical.
The five differences per profile are three path manifests and two Cytoscape.js
producer labels. Explicit normalization of analysis-root paths and the producer
name makes every compared file identical. Scientific attributes are not
normalized. PDF metadata and other summary/status files are not covered by
this report. This is not yet a complete whole-analysis comparison.

A second summary/export generation matched the first byte-for-byte for all
25 v2 and 59 v3 files. The installed-package process confirmed that the
original `hotnetR` namespace was not loaded; this does not substitute for a
clean-machine dependency restoration test.

## Package validation

hotnetR2 0.0.1 passed the complete test suite and metadata-free macOS
`R CMD check --no-manual` with built vignettes: **Status: OK**.
The checked development tarball is
`/Users/li26191/tmp/r-package-hotnetR2.enoIip/hotnetR2_0.0.1.tar.gz`;
the check log is retained alongside it under `hotnetR2.Rcheck/00check.log`.
No CRAN incoming check or PDF manual was requested. No release tag is created.

The latest checked source includes the macOS nonexistent-path symlink fix.
The unneeded copied journal PDF was moved outside the package source before
this build; the original hotnetR asset was left untouched.

An archival bundle of 155 files (about 1.5 GB) was copied and checksum-verified
with tools/export_reference_bundle.R. Its resource manifest is retained here.
It includes native executables/extensions and is not a cross-platform lock or
permission to redistribute the input data. Dependency restore remains pending.

## Independent LDAK and relocation results

Fresh v3 LDAK `remls.all` and `genes.details` match for enhancer and promoter
calculations, and both standardized score tables match after decompression.
Fresh v2 enhancer results and promoter `genes.details` match. However, v2
promoter permutation statistics differ: all observed statistics match, while
17,565 `LRT_P_Perm` values differ (maximum absolute difference 0.000970).
The first per-gene permutation difference is RYR2. This remains an unresolved
external-run reproduction failure, not an accepted tolerance.

Relocating the 155-resource bundle into a project path containing spaces and
running annotation/classification from an unrelated working directory reproduced
all five selected v3 annotation files. The existing R library was used; a
clean-machine software restore has not been demonstrated.

The reusable test scripts, actual commands and extended reports now live in
the sibling `../test_hotnetR2/` folder, as requested by the user.
