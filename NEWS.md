# hotnetR2 0.0.4

- Refresh alias/alias_nodup from the September 2026 curated JEME, HiC-PO and
  HiC-PP workbooks, retaining uncovered legacy keys. Resolve C11orf48 to LBHD1
  and MEGT1 to LY6G6D across every source.
- Ship source workbooks, original Stata/R/Rmd programs, an external portable R
  builder, and provenance/conflict audits. The builder is not run on installation.
- Add create_alias_references() for alias-only cache generation; full reference
  setup now uses the same packaged aliases by default. Existing caches are not
  overwritten unless explicitly requested.
- HiC and ARACNe helpers honor the new explicit raw-to-approved mapping marker;
  unmarked historical caches retain their previous direction. JEME ENSG-first
  harmonization is unchanged.

# hotnetR2 0.0.3

- TSS profiles now retain all three JEME enhancer classes in both LDAK inputs
  and regulatory edges. HiC-PO and gene-body defaults remain class1-only.
- Honor `drop_non_class1_edge_sources` consistently in annotations and networks.
  Class labels and class-specific audit files are retained for all regions.
- Existing explicit `[JEME, HiC_PO]` lists preserve the old class1-only policy.

# hotnetR2 0.0.2

- Empty or omitted regulatory HiC tissue selections now disable HiC in both
  annotations and networks. JEME-only runs do not require HiC caches or aliases.
- Preserve typed empty HiC nodes/edges throughout classification and scoring.
- Existing nonempty v2/v3 HiC settings and get_hic(NULL) behavior are unchanged.

# hotnetR2 0.0.1

- Development migration based on hotnetR 0.4.6.2 and the v2/v3 workflow.
- Prefer v2/v3 decisions over earlier analysis2 behavior.
- Keep legacy public data interfaces; give analysis-specific filtering and
  LDAK summarization distinct names to avoid function collisions.
- Independent package namespace and installed-resource lookup.
- Portable v2/v3 project initialization, checksummed resource provisioning,
  bedtools classification, and installed workflow execution entry points.
- Full-data annotation, score, network and fresh HotNet verification reports;
  clean-machine and fresh LDAK certification remain pending.
- Correct inherited Roxygen/data documentation and external LDAK input paths.
