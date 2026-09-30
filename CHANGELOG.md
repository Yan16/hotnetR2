# Decision-oriented changelog

Consolidated 2026-09-23 against source a722b47. `NEWS.md` remains the release
record. This file explains which decisions remain current; it is not a full
transcript. See PROJECT_SPEC IDs and KNOWN_ISSUES for limitations.

| Version / evidence | Change still in force | Replaced or corrected approach |
| --- | --- | --- |
| 0.2.3, unit tests and real 127-tissue JEME cache expansion | The LDAK summary stage writes the additive `enhancer_ldak_long.tsv.gz`: every standardized enhancer score is retained, JEME enhancers expand to distinct raw `promoterFull`/`ENSG` plus harmonized promoter/tissue/tissue-name rows, and non-JEME enhancers retain one row with missing target fields. The exported function can update an existing summary without rerunning LDAK. | Using the aggregated annotation `tissues` string to infer enhancer targets, dropping HiC-only enhancers, or rerunning association testing merely to create the reporting table. |
| 0.2.2, rendered schema-reference vignette and package check | `analysis.yaml` schema version 1 is documented field by field, separating required/conditional keys, executable fallbacks, canonical-profile examples, legacy aliases and inert/reserved fields. The reference also states which stage must be regenerated after each class of change. No runtime default or scientific setting changed. | Treating every value shown in an example YAML as an implicit default, or assuming accepted but unconsumed fields affect computation. Future removals or changed defaults require a new schema version, validator tests and migration notes. |
| 0.2.1, rendered vignette and example-config checks | The documented canonical analysis is the corrected E092/Gastric TSS workflow with no enhancer-class exclusion and current curated harmonization. v6/v7 supply the four-network and permutation/manual-delta orchestration patterns, but their frozen scientific outputs are not canonical corrected-coordinate baselines. Package-local ARACNe and STRING examples replace copying historical folders. | Treating original analysis2, early v2--v5 folders, gene-window/class-exclusion profiles, or frozen v6/v7 caches as current templates. Historical alternatives remain described briefly for possible explicit future revival. |
| 0.2.0, STRINGdb workflow tests | `interaction_network.type=stringdb` makes Network3/4 STRING-augmented forms of Network1/2. Human STRING v12, inclusive score 700, cached inputs, existing promoter-gene endpoints only, deterministic legacy shortest-name mapping, no STRING loops/duplicate pairs, and complete audits are explicit profile settings. Omitted backend preserves legacy ARACNe behavior. | Reusing ARACNe-specific controls for STRING, adding STRING-only nodes/enhancers, implicit downloads during analysis, or silently accepting a different STRING release. |
| 0.1.1, coordinate regression suite | All active LDAK and bedtools intervals are 0-based half-open. NCBI `+`/`-` TSS and gene-body conversions are strand-aware and explicit; JEME/HiC enhancer intervals remain BED; touching boundaries do not overlap. | The frozen v6/v7 extra start subtraction and positive-strand `START+1` TSS arithmetic. Outputs regenerated with 0.1.1 are corrected scientific results, not byte-identical legacy reproductions. |
| 0.1.0, local coordinate audit | Dataset coordinate/assembly contracts are centralized in `docs/dataset_coordinate_reference.md`. Exact allele-compatible BIM comparisons establish BEEA `POS` and MVP `pos_hg19` as 1-based variant positions. JEME `location` is a 1-based hg19 TSS, JEME enhancer bounds are BED, and active JEME/HiC LDAK promoters are rematched to 1-based-inclusive NCBI gene bounds. | Treating GWAS positions as BED starts, assigning one blanket convention to every JEME/HiC field, or applying an unrecorded conversion. Legacy v6/v7 interval arithmetic is still preserved pending a separately versioned correction; HiC other-end fragment provenance remains unresolved. |
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
