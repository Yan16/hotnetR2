# Reproduction and clean-machine acceptance contract

## Frozen references

| Measure | v2_analysis2 | v3_analysis2 |
|---|---:|---:|
| Enhancer universe | 446623 | 446623 |
| Class 1 | 225571 | 424653 |
| Class 2 | 196458 | 1951 |
| Class 3 | 24594 | 20019 |
| Promoter annotations | 18816 | 18816 |
| Enhancer scores | 212897 | 404183 |
| Promoter scores | 17571 | 17253 |
| Network 1 nodes / edges | 733 / 974 | 907 / 1172 |
| Network 2 nodes / edges | 1170 / 1467 | 2000 / 2528 |
| Network 1 multi-node clusters | 3 | 11 |
| Network 2 multi-node clusters | 2 | 11 |

Counts are smoke checks, not identity proof. The v3 reported permutation
p-value 0.0 means zero extreme permutations among 100 according to the current
HotNet estimator; it does not establish that the underlying probability is zero.

## Resource bundle

The reproduction project includes relative resource paths and checksums for:

- BEEA LDAK summary, p-values and extract (user-supplied GWAS).
- EUR404 LD .bed, .bim and .fam, plus reference provenance.
- Every raw JEME/Hi-C file and tissue-key mapping used by importers.
- NCBI GRCh37 feature table, GENCODE identifier map, and alias map.
- LDAK executable/build, HotNet source/patches, Python lock, bedtools version.
- R package tarball/hash, R version, renv.lock, OS/architecture and locale.

Original caches are reproducibility inputs until regeneration from pinned raw
resources is proven identical. Check licensing/redistribution before placing
assets in package tarballs; manifests may point to user-provided bundles.
No private GWAS data or large LD panels belong in Git or the R package.

A new machine must not require ../analysis2, ../analysis1, an author-specific
/Users path, a preinstalled hotnetR, or implicit here() cache discovery.

## Artifact identity

Capture manifests before code migration. Each entry records analysis, stage,
relative path, size, SHA-256, comparison mode, and exception reason if any.
Gzip files also record decompressed SHA-256 for diagnosis.

Compare all deterministic package data: .loc files, annotation/detail/class
tables, standardized scores, network indexes/edges/scores/metadata, aggregate
and individual cluster tables, and graph exports. Missing and unexpected files
must be failures. Check schemas, order, duplicate multiplicity and values.

Prefer raw-byte equality in a pinned environment. Define exceptions in advance:

- gzip header timestamps: compare decompressed bytes and report raw mismatch;
  make future writers deterministic without rewriting frozen baseline files;
- embedded absolute project paths: normalize only declared path fields, not
  every occurrence of a string or arbitrary table columns;
- logs/timestamps/tool status and PDF metadata: provenance/structure checks;
- version-dependent RDS serialization: exact object equality plus raw hashes.

Record exact exception paths/columns and reasons. Broadly excluding all
summaries or silently accepting numerical tolerance is not permitted.

Two gates are distinct: (1) reproduce package outputs from frozen external
LDAK/HotNet results, (2) reproduce the external calculations themselves.
Passing gate 1 cannot be reported as a complete independent rerun.

## Clean-machine test protocol

1. Install hotnetR2 from the checked tarball into an empty R library.
2. Restore locked dependencies and toolchain; record actual hashes/versions.
3. Provision the checksummed bundle in an arbitrary new root with spaces.
4. Block access to original source trees, libraries, home and caches.
5. Run compact installed-package fixtures, then v2 and v3 in separate outputs.
6. Compare full file inventories and artifact identities, retaining reports.
7. Repeat to verify deterministic outputs and valid resume.
8. Change a fixture input and confirm downstream invalidation.
9. Record environment image digest or clean-machine evidence and commands.

Cross-platform exact numeric identity is not promised solely from a seed:
LDAK builds, CPU architecture, BLAS, Python dependencies, compression and
serialization can differ. Certify specific environments. If original binary
reproduction is impossible, document the limitation and keep native reruns
separate from the exact compatibility certificate.
