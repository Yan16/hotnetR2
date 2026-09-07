# hotnetR2 implementation plan

## Objective and scope

Create an independently installable hotnetR2 package from hotnetR 0.4.6.2 and
the v2/v3 analysis helper implementations. Reproduce both complete reference
analyses from a portable input bundle on a machine without this workspace.
Preserve useful existing API signatures, consolidate duplicate decisions,
prefer explicit tidyverse calls, and document exports with Roxygen.

Initial release scope is v2 and v3 regulatory Networks 1–2. Existing hotnetR
APIs needed by these workflows are ported; remaining public APIs are inventoried
and either preserved or explicitly deferred. Network 3/4 and additional cohorts
are future coverage, not implied release guarantees.

## Phase 0 — Freeze evidence and resolve contracts

Capture source commits, working-tree hashes, installed package versions,
reference configs, data files, tool versions, and output manifests. Commit
the workflow, architecture decisions, API map, and detailed checklist first.

Investigate the current interval inconsistency: NCBI parser retains input
starts; v3 adds one on the plus strand under a BED assumption; classifier
subtracts another base before bedtools; LDAK reads .loc as BED. Freeze actual
v2/v3 arithmetic as named compatibility profiles. A corrected biological
coordinate policy needs independent boundary tests and a separately named run.

Gate: all inputs and deterministic outputs have explicit identities, all
API collisions have destinations, and interval policies are documented.

## Phase 1 — Standalone package foundation

Create DESCRIPTION, license attribution, Roxygen namespace, focused modules,
tests, installed scripts, and vignettes. Port the dependency closure from
hotnetR and helper code, keeping the source repository unmodified.
Remove calls to hotnetR:: and hotnetR::: and analysis-local source() bootstraps.
Preserve data schemas, row order, duplicate behavior, and output formatting.

Gate: installed hotnetR2 loads with original hotnetR unavailable; migrated
API formals are compared against recorded signatures.

## Phase 2 — Portable inputs and environment

Create project-relative resource manifests with resource IDs, SHA-256, version,
assembly, coordinate convention, URL or user-supplied status, and producer.
Bundle small lawful mapping assets and parser fixtures; transfer large/private
inputs separately with checksums. Pin JEME tissue keys, Hi-C supplementary data,
NCBI, GENCODE, and alias mappings. Avoid reconstructing old mappings from
unversioned current services.

Pin R dependencies in a reproduction-project renv.lock and Python dependencies
in a locked environment. Record HotNet source commit and compatibility shim,
LDAK binary hash/build/architecture, and bedtools version. Replace personal
absolute paths and ambient here()/working-directory cache discovery.

Gate: restore a compact project at a path containing spaces, in a clean
library and empty home/cache, using installed resources only.

## Phase 3 — Annotation and classification

Unify all-tissue JEME ENSG mapping and Hi-C symbol mapping; use one saved map
for annotations and tissue-specific edges. Implement explicit gene-body and
TSS profiles with frozen legacy coordinate semantics. Bring the Bash/bedtools
classifier behind an R function; use pure R/tidyverse for tables and bedtools
for large interval intersections. Use portable hashing/sorting and safe quoting.

Remove runtime dependence on analysis2 source files: generate the same complete
enhancer universe from pinned resources. Separate source-edge comparison
audits from the scientific classification stage.

Gate: exact .loc, classification, harmonization and details-table comparisons
against both references; boundary/strand/containment fixtures pass.

## Phase 4 — LDAK and scores

Unify buffer handling and p-value transformation. Preserve legacy
summarize_ldak_results() arguments; give config-driven summarization an explicit
new name. Execute checked commands with captured status and log paths.
Check result products in addition to exit status because LDAK has reported
input errors with successful shell status in the existing workflow.

Gate: regenerated summaries match existing raw results; new isolated LDAK
runs are separately compared under pinned executable/environment settings.

## Phase 5 — Networks, HotNet, and exports

Reuse canonical harmonization/classification output when constructing E094
and gastric Hi-C edges. Preserve self-loops, row ordering, duplicate policy,
node indexing, component filtering, and score caps. Keep ARACNe disabled for
these two profiles. Port cluster parsing, per-cluster integrity, and exports.

Replace file-presence resume with stage fingerprints of inputs, config, code,
and tools. Invalidate downstream products when fingerprints change.

Gate: exact network TSVs, clusters, and deterministic graph exports; complete
observed plus 100 permutation products; repeated execution reproduces data.

## Phase 6 — Clean-machine reproduction and release

Run compact tests from an installed tarball, then both full profiles in fresh
output roots. Restore and execute on a separate clean VM/container or machine.
Provide an archival route for exact historical resources/executables and a
native route for supported platforms with explicitly limited numerical claims.

External floating-point behavior can differ across CPU/BLAS/tool builds.
Require exact outputs in the pinned compatibility environment; report native
cross-platform differences rather than silently adding tolerance or accepting
matching counts alone. Public data may be downloadable; controlled GWAS must
be supplied by the authorized user.

Gate: full tests, built vignettes, R CMD check Status: OK, installed-package
tests, complete artifact reports, clean-machine evidence, and documented
restore/run commands. Tag v0.1.0 only after all gates pass.

## Coding and review strategy

Keep import validation at public/file/tool boundaries. Internal package-owned
tables follow documented contracts; do not repeatedly validate every column
or catch errors and substitute empty results. One transformation per function,
explicit joins, explicit relationship/cardinality where relevant, and stable
ordering at serialization boundaries. Use dplyr::, tidyr::, readr::, purrr::,
tibble::, stringr:: etc.; base and stats calls may remain unqualified.

Each phase is a focused Git commit with tests/evidence and updated task status.
Never commit input cohorts, installed libraries, caches, or generated large
analyses. Separate compatibility refactors from scientific corrections.
