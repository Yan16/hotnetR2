# External alias-data preparation

This directory is shipped for reproducible data preparation, not sourced by
the package on load or installation. The original `.do`, `.R` and `.Rmd`
programs are preserved byte-for-byte in
`extdata/harmonization_2026_0909/original/`. Source workbooks are adjacent.

From an installed package:

```r
source(system.file("scripts", "harmonization", "build_alias_data.R", package="hotnetR2"))
build_alias_data(
  source_dir = system.file("extdata", "harmonization_2026_0909", package="hotnetR2"),
  legacy_file = system.file("extdata", "gtex_hugo_merged_by_ens_SORTens_FINALv1.xlsx", package="hotnetR2"),
  output_dir = "rebuilt_aliases"
)
```

Requires readxl, dplyr, tibble, stringr, readr, purrr and digest (package
dependencies). The builder writes RDA/RDS datasets and TSV audits. It uses
the delivered final workbooks as authoritative inputs, not live HGNC services.
It preserves the legacy 13-column schema and uncovered keys; curated keys
replace legacy target symbols/metadata. New keys have NA for unavailable
legacy/ENSG fields. Unavailable curated Entrez IDs are not copied from stale
legacy target metadata.

HiC-PP contributes both endpoint lookup keys, not repeated network edges.
Conflicts are checked across all sources after applying the two project
decisions: C11orf48 -> LBHD1 and MEGT1 -> LY6G6D. No rule remaps C11orf98 itself,
and no rule globally replaces LY6G6F-LY6G6D. Further unresolved conflicts stop
the build. Classifications and ENSG harmonization are outside this builder.

The legacy Stata scripts require their original external inputs, Stata and
Windows/network paths. The original R scripts are retained as historical
preprocessing references, not certified equivalents of those Stata programs:
for example, HiC's chromosome-sensitive Excel-date repairs in Stata are more
specific than the original R lookup. The final delivered workbooks already
encode that curation. No upstream Stata rerun is required to rebuild these
package datasets. The original Rmd reports raw workbook conflicts; consult
`manual_resolution_audit.tsv` for the subsequent user decisions.

The generated tables carry `mapping_direction = "raw_to_approved"` and a
release attribute. Do not strip those attributes when saving RDS caches.
The raw source files are unchanged; `source_manifest.tsv` records input hashes.
