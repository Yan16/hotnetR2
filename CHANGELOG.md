# Decision-oriented changelog

Consolidated 2026-09-23 against source a722b47. `NEWS.md` remains the release
record. This file explains which decisions remain current; it is not a full
transcript. See PROJECT_SPEC IDs and KNOWN_ISSUES for limitations.

| Version / evidence | Change still in force | Replaced or corrected approach |
| --- | --- | --- |
| 0.0.1, e64783f; earlier cf1353e/068bc48 | Independent hotnetR2 namespace and sibling repository; installed-resource lookup; explicit profile configuration; resource manifests; source-specific names avoid function collisions. | Sourcing arbitrary analysis-local hotnethelper copies as runtime dependencies. |
| v2/v3 migration, documented tests | Enhancer expansion 1 kb; promoter expansion 2 kb; three enhancer classes; all-tissue ENSG-first JEME mapping; optional HiC/ARACNe; gene-body versus TSS profiles. | Zero-flank enhancer classification, two-class overlap model, JEME symbol-only mapping. Gene-body remains supported, not globally removed. |
| 0.0.1 verification | External LDAK paths handled, package data/Roxygen corrected, source assets located from installed package, result writes isolated. | Implicit original repository and author-library assumptions. Clean-machine certification is still incomplete. |
| Current cluster-table helper; test-hhotnet-results.R | Per-cluster outputs subset by the actual requested cluster ID; duplicate node membership is rejected. | Historical per-cluster files containing multiple IDs. Original fixes also affected old analysis folders; do not infer all old exports were regenerated. |
| 0.0.2, ab2f051 | Empty/absent workflow HiC tissue disables HiC; typed empty tables supported through annotations/networks. | Mandatory HiC selector and cache use in JEME-only workflows. Public get_hic(NULL) retains all-tissue semantics. |
| 0.0.3, 3b75a7a | TSS defaults retain every JEME enhancer class; explicit drop list governs both annotation and edges. | TSS universal class1-only policy. HiC PO remains class1-only by default. |
| 0.0.4, d31b460, fe901f7, 23afc28 | September curated workbooks, explicit two alias resolutions, mapping-direction attributes, external builder and original program-byte preservation. | Incomplete packaged alias tables; ambiguous map direction for new data. Legacy cache behavior is deliberately preserved. |
| 0.0.5, a56405f | HiC edge-type selection before annotation/network construction; PP-only mode needs no PO cache. | Dropping PO only in final graphs while retaining exclusive PO annotations. |
| 0.0.6, cd43f18 | Optional HiC-PO exclusion against selected JEME tissue union; both flank defaults 1 kb; separate overlap audits. | Comparing against all JEME tissues indiscriminately, or conflating the new filter with TSS classes. |
| 0.0.7, a722b47 | Similarity-only HHN stage; observed-only hierarchy, vectorized manual delta interface, native cut semantics, delta-named TSV/CX2, no permutation significance. | Requiring null hierarchies to obtain exploratory cuts. Standard permutation mode still supported. |
| 0.0.7 packaging | Exclude generated compat Python bytecode from source tarballs; metadata-free macOS build/check. | Shipping machine-generated Python cache as release source. Tracked historical bytecode cleanup is still pending. |

Project-only launcher changes (minimum-version acceptance, folder-derived names,
user-library loading) are documented in the sibling project's changelog, not
claimed as package API changes. No global promoter-window replacement, cache
refresh, or new statistical method is authorized by this consolidation.

## Verification status

- 0.0.7 release: full tests, R CMD check with vignettes (`Status: OK`), temporary
  installed R4.5 tests using real local HHN on synthetic inputs. Saved in
  `../test_hotnetR2/reports/MANUAL_DELTA_0.0.7.md` and the corresponding check log.
- Historical v2/v3 selected artifact identity and later LDAK variability are
  recorded in `docs/verification.md`; no blanket exact-reproduction claim.
- 2026-09-23 consolidation: read-only v6/v7 code/config/result audit and new
  checksum evidence; no new scientific computation, release or installation.
