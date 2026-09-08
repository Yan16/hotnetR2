# Detailed implementation tasks

Status legend: checked means completed with evidence; unchecked means pending.
A passing later gate does not waive unfinished earlier tasks.

Implementation update (2026-09-08): v2/v3 take precedence over analysis2.
See docs/verification.md for the exact scope of full-data comparisons; unchecked
full-folder/clean-machine gates remain required even where selected files match.

## P0 — Design and frozen baseline

- [x] P0.01 Inspect hotnetR metadata, public functions, v2/v3 configuration and scripts.
- [x] P0.02 Create workflow schematic, implementation plan and module ownership.
- [x] P0.03 Record API collisions, coordinate inconsistency and resume limitation.
- [x] P0.04 Define coding, Roxygen and clean-machine acceptance contracts.
- [x] P0.05 Initialize package Git history and commit design documents.
- [ ] P0.06 Record function signatures/exports and migration disposition for every API.
- [x] P0.06a Capture 323 function definitions/exports across 50 source/metadata files;
  report duplicate definitions and conflicting signatures automatically.
- [ ] P0.06b Assign and review a migration disposition for every captured public API.
- [ ] P0.07 Capture working source SHA-256, installed package versions and session info.
- [ ] P0.08 Hash frozen inputs and reference outputs; classify exact comparison exceptions.
- [ ] P0.09 Audit source NCBI coordinates/assembly records and quantify boundary differences.
- [ ] P0.10 Freeze v2/v3 legacy profiles with expected outputs and source manifest IDs.
- [ ] P0.11 Document asset licenses, distribution eligibility and user-supplied inputs.
Gate: reviewed baseline manifest, complete API map, explicit coordinate policies.

## P1 — Package foundation and API migration

- [x] P1.01 Create DESCRIPTION, MIT attribution, NEWS, .Rbuildignore and namespace.
- [x] P1.02 Port required hotnetR importers, mappings, graph and score functions.
- [x] P1.03 Port reconciled helper config, annotation, LDAK and network implementations.
- [x] P1.04 Resolve duplicate exclude_MHC, summarize_ldak_results and public/internal
  filter_network_ldak interfaces explicitly.
- [x] P1.05 Replace hotnetR::/::: and hardcoded package lookup with internal calls.
- [x] P1.06 Remove analysis-local source() bootstraps from installed CLI entry points.
- [ ] P1.07 Apply explicit tidyverse namespaces; keep internal contracts concise.
- [ ] P1.08 Add Roxygen params, return schemas, coordinate rules, examples and side effects.
- [ ] P1.09 Generate NAMESPACE/man; compare compatible public formals automatically.
- [ ] P1.10 Install and load without original hotnetR or any analysis source directory.
Gate: coherent installed package, no API collisions or hidden runtime source dependency.

## P2 — Resource provisioning and portable projects

- [ ] P2.01 Define resource manifest schema and checksummed resource IDs.
- [x] P2.02 Add initialize_analysis() templates for v2/v3 with relative paths.
- [x] P2.03 Add explicit cache/resource context to readers without breaking leading args.
- [ ] P2.04 Freeze tissue keys, NCBI, GENCODE and aliases with provenance.
- [ ] P2.05 Implement verified URL downloads and user-supplied input registration.
- [ ] P2.06 Support exporting/importing archival resource bundles outside the tarball.
- [ ] P2.07 Reproduce caches from raw files and compare to frozen cache objects.
- [ ] P2.08 Record LDAK executable build/hash and bedtools version.
- [ ] P2.09 Pin HotNet checkout, Python dependencies and compatibility shim hashes.
- [ ] P2.10 Supply reproduction-project R lockfile; test clean restore.
- [ ] P2.11 Verify paths with spaces and no implicit home/cache dependencies.
Gate: minimal project restoration using only installed package and manifest resources.

## P3 — Harmonization, intervals and enhancer classes

- [ ] P3.01 Implement one ENSG normalizer and all-tissue mapping engine.
- [ ] P3.02 Preserve ambiguous/unmapped fallbacks and Hi-C mapping semantics.
- [ ] P3.03 Preserve multi-promoter splitting, tissue membership and stable ordering.
- [x] P3.04 Implement explicit legacy gene-body and TSS interval adapters.
- [ ] P3.05 Implement canonical BED conversion separately with source-verified tests.
- [x] P3.06 Move overlap/class precedence into one classify_enhancers() operation.
- [x] P3.07 Package bedtools invocation with portable sort/hash tools and quoting.
- [x] P3.08 Generate source enhancer universe without reading analysis2 outputs.
- [ ] P3.09 Separate legacy source-edge audits from runtime classification outputs.
- [ ] P3.10 Test no overlap, endpoint touch, full containment, partial overlap, multiple
  windows, both strands, chromosome mismatch and chromosome-start clipping.
- [ ] P3.11 Compare complete annotation/classification artifacts for both profiles.
Gate: matching .loc/details/classes and explicit audit of any policy differences.

## P4 — LDAK execution and summaries

- [x] P4.01 Generate cut/calc/join commands from canonical config.
- [ ] P4.02 Centralize buffers, selected SNP inputs and executable resolution.
- [ ] P4.03 Validate tool outputs as well as exit status; preserve error logs.
- [ ] P4.04 Preserve legacy summarize_ldak_results() formals and behavior.
- [x] P4.05 Implement summarize_analysis_ldak_results() and update callers.
- [ ] P4.06 Centralize finite/infinite score transformation and FDR handling.
- [x] P4.07 Compare summaries from existing frozen LDAK raw products.
- [ ] P4.08 Rerun both LDAK profiles independently and compare full outputs.
Gate: structural integrity and exact package-output identity, separate external-run report.

## P5 — Networks and HotNet execution

- [ ] P5.01 Use saved canonical mapping and enhancer classes in network construction.
- [ ] P5.02 Preserve node order, row multiplicity, self-loops and component filters.
- [x] P5.03 Keep disabled ARACNe explicit and restrict current profiles to Networks 1–2.
- [x] P5.04 Compare full nodes/scores/index/edge/metadata tables against both references.
- [ ] P5.05 Implement stage fingerprints and dependency-aware invalidation.
- [ ] P5.06 Package local Python execution and portable progress/status reporting.
- [x] P5.07 Run observed + 100 permutations with recorded environment.
- [ ] P5.08 Verify clean resume and invalidation after a real fixture input change.
Gate: network artifact identity and complete validated HotNet runs.

## P6 — Cluster outputs and artifact verification

- [x] P6.01 Port corrected cluster parsing and per-cluster splitting.
- [x] P6.02 Port deterministic TSV/GraphML/Cytoscape.js/CX2 serializers.
- [ ] P6.03 Implement missing/unexpected/schema/order/byte comparison reports.
- [ ] P6.04 Apply only predeclared path/compression/runtime metadata exceptions.
- [x] P6.05 Verify cluster IDs and all exported edge endpoints.
- [x] P6.06 Regenerate summaries and exports twice and confirm deterministic identity
  for the declared HotNet comparison inventory (25 v2 / 59 v3 files).
- [ ] P6.07 Keep raw p=0 reporting with permutation-count interpretation.
Gate: complete artifact inventory and zero unexplained differences.

## P7 — Documentation, independent validation and release

- [ ] P7.01 Provide installed-package quickstart and full v2/v3 run guides.
- [ ] P7.02 Build workflow, migration, resource and reproduction vignettes.
- [x] P7.03 Run focused tests, full suite and macOS metadata-free R CMD check
  (development 0.0.1, Status: OK; see docs/verification.md).
- [ ] P7.04 Install checked tarball into a clean library and run installed tests.
- [ ] P7.05 Restore compact fixtures on a clean VM/container or separate machine.
- [ ] P7.06 Execute both full analyses without original paths/caches/source trees.
- [ ] P7.07 Produce per-platform artifact reports with toolchain/environment identities.
- [ ] P7.08 Verify second-run identity and input-change invalidation.
- [ ] P7.09 Publish exact support/limitations and user-supplied resource requirements.
- [ ] P7.10 Update project README index and task statuses with concrete evidence.
- [ ] P7.11 Commit release notes and tag v0.1.0 only after all required gates pass.
Gate: independently reproducible analyses, checked installed release and retained evidence.

## Suggested commit sequence

1. docs: define workflow, compatibility contracts and implementation roadmap
2. test: freeze source/API/resource and artifact baselines
3. feat: add standalone package and compatible core API
4. feat: add portable resource manifests and project initialization
5. refactor: unify harmonization and explicit interval policies
6. feat: add enhancer classification and LDAK execution
7. feat: migrate regulatory networks and fingerprinted HotNet runner
8. test: certify cluster/export identity and clean-machine restoration
9. docs: finalize reproducible release instructions and evidence
