# Source inventory

Generated from the existing hotnetR and the v2/v3 helper source trees before
implementation. Run from the hotnetR2 repository root:

```bash
R_LIBS_USER="$PWD/../buas_paper_2026/v3_analysis2/Rlib" \
  Rscript tools/capture_source_inventory.R "$PWD/../buas_paper_2026" "$PWD/../hotnetR"
```

The script needs dplyr, purrr, readr, stringr, tibble and digest. The existing
library above is a development convenience for this inventory only; the
finished package must support an independent clean library.

| Report | Meaning |
|---|---|
| source_files.tsv | 50 R source, DESCRIPTION and NAMESPACE files with raw SHA-256. |
| function_signatures.tsv | 323 top-level definitions, full formal arguments, export status and body hashes. |
| duplicate_definitions.tsv | Same-source duplicate function definitions. |
| signature_conflicts.tsv | Names appearing across sources with differing formal arguments. |

No source is executed or loaded. Function-body hashes use serialized R
expressions and are diagnostic; raw file hashes identify the source bytes.
Reports deliberately exclude absolute machine paths and capture timestamps
so a repeated inventory on unchanged inputs is deterministic.

This is a source baseline. It does not replace the pending raw-input and
analysis-output manifests, or prove that the installed hotnetR binary matches
the current source tree.
