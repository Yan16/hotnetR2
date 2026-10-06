# Known issues and refactoring gates

As of 2026-09-23, runtime baseline a722b47 / 0.0.7. “Confirmed” below means
observed code/artifact behavior or previously recorded reproducible failure;
it does not automatically authorize changing scientific outputs.

| ID / priority | Status and evidence | Required next step / acceptance |
| --- | --- | --- |
| K01 / resolved in 0.1.1 | **Resolved mixed-coordinate conflict.** NCBI and JEME TSS values are converted once to BED; NCBI gene bodies are converted once; JEME and Jung HiC enhancer/restriction-fragment intervals remain BED. Classifier and HiC-JEME overlap helpers no longer subtract an extra base. | Regression tests cover `+`/`-` TSS, NCBI gene body, JEME/NCBI public converters, zero starts, touching boundaries and BED-based MHC overlap. Results made with 0.1.1 are intentionally coordinate-corrected and are not byte-identical legacy v6/v7 reproductions. |
| K02 / high | **Confirmed historical failure:** fresh v2 promoter LDAK permutation p-values differ even with same seed/binary/8 threads; observed statistics match. `docs/verification.md`: 17,565 reference differences, another 17,570 between repeat runs. | Diagnose external threading/RNG; no invented tolerance. This does not prove v6/v7 will fail, nor prove their fresh reproduction. |
| K03 / high | **Confirmed baseline/cache divergence:** active v6/v7 alias cache is unmarked legacy data, whereas fresh package setup defaults to marked September-curated tables. | Restore frozen cache for exact target reproduction. Separately test curated refresh as a scientific change; never silently overwrite it during setup. |
| K04 / high | **Not certified:** full dependency restoration and complete fresh v6/v7 reproduction on a new computer. Downloads may change, native tools differ, and permissions/licenses constrain redistribution. | Freeze source inputs, R/Python locks, native compilation and exact versions; execute isolated restoration before declaring portability. |
| K05 / medium | **Confirmed:** initialize_analysis supports only v2/v3 templates; generic defaults differ from target v6/v7 overrides. Project launcher and package config decisions are split. | Add explicit schema/profile resolution with tests; preserve old API names and historical profiles. |
| K06 / medium | **Code-review risk, not demonstrated failure:** invalidation/resume is heterogeneous. Manual fingerprints omit some runtime environment identity compared with permutation execution; direct extraction accepts existing matrix files without proving they match current network inputs. | Characterization tests changing inputs, Python/native libraries, overrides and stale matrices. Add dependency-aware stage signatures before claiming robust resume. |
| K07 / medium | **Confirmed:** manual-grid manifest describes latest invocation only; old delta files persist. Numeric filename labels can collide (explicitly rejected); clusters can be singleton-only. | Define run/grid identity and inventory policy; preserve empty-graph handling and do not delete old outputs implicitly. |
| K08 / medium | **Confirmed:** standard CX2 selection includes surrounding connected components; manual CX2 contains only qualifying cluster members and within-cluster edges. | Keep separate documented strategies. Never refactor both through one default that changes node membership. |
| K09 / medium | **Confirmed source behavior:** local permutation seed includes base+i; Apptainer branch currently uses i directly. Slurm HHN script hardcodes 128G despite project slurm.memory=16G. | Explicitly test/fix only with scope approval; local v6/v7 do not use Apptainer/Slurm submission. Clarify effective resource semantics. |
| K10 / medium | **Confirmed limitation:** analysis annotation universe supports ALL; selected-tissue network/filter references are separate. Annotation presence is not network participation. | Keep distinction and add explicit universe types; don't advertise tissue-restricted annotation execution as supported. |
| K11 / medium | **Confirmed provenance gap:** final alias workbooks can be rebuilt into package tables, but original Stata and original R curation are not certified equivalent; Stata depends on external paths/data. | Preserve originals and provenance. Test portable builder against final workbooks, not an assumed upstream rerun. |
| K12 / low | **Confirmed:** old Python `.pyc` remains tracked in source history, but excluded from 0.0.7 tarball; setup script messages still mention hotnetR. | Hygiene/documentation cleanup in separate commits; preserve originals where byte identity is required. |
| K13 / medium | **Unverified edge cases:** normalized punctuation keys may collapse distinct IDs; very small/degenerate graphs, missing finite scores, new tissue/chromosome variants and same-name/different-coordinate intervals need broader characterization. | Add external-input fixtures before simplifying validation. Distinguish expected singleton cuts from empty/invalid hierarchy input. |
| K14 / low | **Confirmed API presentation gap:** validate_hhotnet_completion in manual mode certifies matrices only, not extraction; some Rd descriptions still emphasize permutation results. | Split or clarify stage-specific completion contracts without claiming matrix presence proves successful cuts. |
| K15 / medium | **STRING identity-mapping limitation:** STRING protein aliases can map several graph gene names to one protein. The 0.2.0 backend deliberately reproduces the historical shortest-gene-name rule with a lexical tie-break and writes conflicts, but this is a deterministic compatibility policy rather than a biological disambiguation. Results also depend on the exact cached STRING release. | Preserve the v12 cache and audit tables with each run. Review mapping conflicts before interpretation; changing the mapping policy or release is a new analysis profile, not a transparent cache refresh. |
| K16 / medium | **Historical profile boundary:** frozen v6/v7 provide validated orchestration evidence but include older classification/cache/coordinate behavior. Original analysis2 and early v2--v5 folders contain superseded combinations. | Use the 0.2.1 package examples for new corrected TSS/no-class analyses. If gene-window or enhancer-class exclusion is revived, give it a separate profile name, tests and expected results; never copy an obsolete folder as the current template. |
| K17 / resolved in 0.2.4 | **Resolved ambiguous summary identifier and lost source-label issue.** Standard LDAK summaries formerly wrote `gene` as an exact duplicate of LDAK's `Gene_Name`, while enhancer outputs retained only harmonized promoter symbols and could not report the original JEME promoter label. | `Gene_Name` is now the sole standardized-summary node identifier and all current summary consumers read it. Both enhancer summary forms contain `original_promoter`: sorted semicolon-collapsed labels in the one-row-per-enhancer table and one raw label per association in the long table. HiC-only rows use `NA`; legacy conflicting aliases fail rather than being silently discarded. |
| K18 / resolved in 0.2.5 | **Resolved missing association-confidence issue.** `enhancer_ldak_long.tsv.gz` formerly omitted JEME `conf_score`, forcing reports to reopen source caches to recover strongest-target scores. | Long summaries now retain source `conf_score` per enhancer--target--tissue association; HiC-only rows use `NA`. The standard enhancer summary intentionally remains unchanged because one enhancer can have multiple target-specific scores. |

## Refactoring order (after fresh-session review)

1. Freeze/audit resources and baseline tests; retain current scientific behavior.
2. Centralize explicit configuration defaults/profile resolution; define interval
   representations and source-specific retention as separate policies.
3. Separate data acquisition, curation, harmonization, annotation, score joining
   and graph construction. Keep setup opt-in and caches immutable by default.
4. Share graph serialization and cluster membership checks while preserving the
   standard/manual export selection distinction.
5. Unify stage manifests, version/environment signatures and resumability.
6. Add installation/resource-bundle and clean-machine tests. Only then simplify
   compatibility shims; retain public names or provide explicit deprecations.

Do not “fix” K02 or K03 by changing inputs/tolerances inside an
otherwise structural refactor. A clean R CMD check is necessary, not proof of
scientific or whole-project reproducibility.
