# Analysis workflow

## End-to-end schematic

```mermaid
flowchart TD
  A[Portable project config + resource manifest + environment locks] --> B[Validate external resources once]
  B --> C[All 127 JEME tissues: normalize ENSG and build one map]
  B --> D[All Hi-C tissues: symbol and alias harmonization]
  B --> E[NCBI GRCh37 gene coordinates and strand]
  C --> F[Harmonized target-gene universe]
  D --> F
  F --> G[Resolve gene coordinates under explicit compatibility policy]
  E --> G
  G --> H{Promoter profile}
  H --> I[v2: full gene interval]
  H --> J[v3: strand-aware TSS interval]
  C --> K[All-tissue enhancer universe]
  D --> K
  I --> L[Apply enhancer 1000 bp and promoter 2000 bp flanks]
  J --> L
  K --> L
  L --> M[Same-chromosome intersection and containment]
  M --> N[Class 1: no overlap]
  M --> O[Class 2: contained in at least one window]
  M --> P[Class 3: remaining partial overlaps]
  N --> Q[Class-1 enhancer LDAK annotations]
  I --> R[Promoter LDAK annotations]
  J --> R
  B --> S[BEEA summary + p-values + extract + EUR404 LD]
  Q --> T[LDAK cut-genes / calc-genes-reml / join-genes-reml]
  R --> T
  S --> T
  T --> U[Standardized enhancer and promoter scores]
  C --> V[E094 JEME using saved all-tissue map]
  D --> W[Gastric Hi-C and split promoter endpoints]
  V --> X[Drop class 2/3 JEME and Hi-C-PO edges]
  W --> X
  N --> X
  U --> Y[Match scores and preserve defined node/edge order]
  X --> Y
  Y --> Z1[Network 1: component max score >= 4]
  Y --> Z2[Network 2: component max score >= 3]
  Z1 --> AA[Hierarchical HotNet: observed + 100 permutations]
  Z2 --> AA
  AA --> AB[Cluster summaries / individual cluster tables]
  AB --> AC[GraphML / Cytoscape JSON / CX2 / TSV]
  AC --> AD[Artifact identity and structural verification]
```

Class 2 has precedence when one window fully contains an enhancer and another
partially overlaps it. Overlap is tested against all same-chromosome target
windows, independently of the enhancer's linked JEME target. Hi-C-PP edges
remain eligible. Both profiles disable ARACNe.

## Module ownership

| Scientific decision | Single owner | Consumers |
|---|---|---|
| Resource identity and location | resources | importers and provenance |
| Identifier policy and ambiguity | harmonization | annotations and network edges |
| Interval meaning and flanks | intervals | classifier, LDAK writer, audits |
| Enhancer class precedence | classification | annotation filter and edge gate |
| Regional p-value to score | scores | node construction and exports |
| Ordering and duplicate policy | network construction | numeric indexes and graphs |
| Resume validity | execution | LDAK and HotNet stages |
| Equality policy | reproducibility | regression tests and release reports |

The interval module must explicitly represent legacy differences between
classification coordinates and LDAK interpretation; a diagram of the desired
flow does not establish that historical coordinate arithmetic was consistent.
See decisions.md before implementing interval conversion.

## Runtime sequence

`initialize project -> provision inputs/tools -> import -> harmonize ->
annotate -> classify -> LDAK -> summarize -> networks -> HotNet -> export ->
verify`.

The finished installed package provides the sequence through documented public
functions and installed command-line entry points. Pure transformations return
tables; orchestration handles files, tools, stage fingerprints, and provenance.
