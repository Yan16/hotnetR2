# Dataset and coordinate-system reference

Status: current coding reference for hotnetR2 0.1.1, audited 2026-09-24.

This document records the datasets that directly define the current BEEA
v6/v7 analysis profiles, plus the separately inspected MVP GWAS input. It is a
coding-boundary reference, not permission to redistribute private data and not
a claim that every upstream raw-data generation step is reproducible.

## Coordinate vocabulary

- **1-based variant position:** a single genomic base numbered from 1. GWAS
  `CHR/POS`, `chr:position` predictor IDs and PLINK BIM base-pair positions use
  this convention. A position must not be decremented when forming an ID.
- **1-based inclusive interval:** both `START` and `END` name included bases.
  An interval `[s,e]` has width `e-s+1`.
- **0-based half-open interval (BED):** `START` is included and `END` is
  excluded. The 1-based inclusive interval `[s,e]` becomes `[s-1,e)`.

LDAK requires the first four `.loc`/gene-file fields to contain name,
chromosome, 0-based start and half-open end. Genomic source intervals must be
converted exactly once at this boundary. GWAS variant positions are identifiers,
not interval starts, and do not receive this conversion.

Primary format references:

- [PLINK `.bim` format](https://www.cog-genomics.org/plink/2.0/formats#bim):
  field 4 is explicitly a 1-based base-pair coordinate.
- [LDAK gene-file format](https://dougspeed.com/ldak-gbat/): gene annotations
  are explicitly 0-start, half-open.
- [NCBI GFF3 coordinate documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/):
  feature start and end are 1-based coordinates. The cached NCBI feature-table
  values are retained as supplied and follow the corresponding feature bounds.
- [Ensembl coordinate documentation](https://grch37.ensembl.org/info/docs/api/general_instructions.html):
  Ensembl coordinates are 1-based and inclusive. This contract applies to the
  Ensembl-derived TSS embedded in JEME `promoterFull` identifiers.

## GWAS coordinate audit

### BEEA_Schroder

Files inspected:

- `GWAS/BEEA_Schroder/head_Meta_beea_20201126_miss2.tsv`, a ten-record subset
  of the unavailable-on-this-machine original
  `Meta_beea_20201126_miss2.tsv.gz`.
- `GWAS/BEEA_Schroder/ldak_BEEA_META.tsv`, the processed LDAK summary input.
- `LD_ref/g1k_EUR404_LDAK.bim`, the GRCh37 reference used by the analysis.

The subset supplies `CHR`, `POS`, alleles, rsID and a `MarkerName` containing
the same coordinate. Eight of ten records have an exact, allele-compatible
`CHR:POS` match in the BIM; the other two are absent from that reference. There
are no position or allele-compatible matches at `POS-1` or `POS+1` for these
records. Processed predictor IDs retain the source value verbatim, for example
source `CHR=10, POS=100000625` becomes `Predictor=10:100000625` and matches the
BIM at base-pair position 100000625.

Conclusion: BEEA `POS`, the coordinate in `MarkerName`, processed `Predictor`,
and the matching BIM base-pair coordinate are **GRCh37 1-based variant
positions**. Do not subtract one. Allele orientation is a separate operation.

Audit identities:

| File | SHA-256 |
| --- | --- |
| BEEA ten-record source subset | `5e60126f687ea757bdf95634fb4386ae49debe9548fa8fa8bdd2be4f15398bdb` |
| BEEA processed LDAK summary | `72803eaf32fee16b53e29d95af4237e42e7ef9600de6a0a077845faadc5c731b` |

The processed summary has 6,134,968 data rows. The accompanying extraction
list has 6,086,498 predictors. These counts and exact files are reproduction
inputs; the incomplete raw-to-processed BEEA provenance remains a later project
TODO and does not change the coordinate conclusion.

### MVP BE-only

File inspected:
`GWAS/MVP/FullSummary_GWAS_BEonly_MVP_hg19_MAF_Higher_01pct.tsv.gz`.

The table contains both `pos_hg38` and lifted `pos_hg19`. The active MVP LDAK
workflow uses `chromosome` plus `pos_hg19`, not `pos_hg38`, and constructs
`Predictor` as `chromosome:pos_hg19` without coordinate conversion.
The retained generation notebook explicitly converts the original hg38
`base_pair_location` to a one-base BED interval `[position-1,position)` before
UCSC liftOver, then assigns the lifted BED end to `pos_hg19` to restore a
1-based position. It renames the unchanged original `base_pair_location` to
`pos_hg38`. Thus both position columns are intentionally 1-based.

For the first 100,000 records, comparison with the GRCh37 BIM found:

| Candidate coordinate | Position matches | Allele-compatible matches |
| --- | ---: | ---: |
| `pos_hg19 - 1` | 1,404 | 288 |
| `pos_hg19` | 88,152 | 88,116 |
| `pos_hg19 + 1` | 1,511 | 291 |

The few offset matches are neighboring reference variants; exact-coordinate,
allele-compatible matching overwhelmingly identifies the convention. The full
existing preprocessing audit independently reports 9,130,225 input records,
6,827,037 exact BIM coordinate matches after QC, and 6,827,033 allele-consistent
final variants.

Conclusion: MVP `pos_hg19` is a **GRCh37/hg19 1-based variant position**.
`pos_hg38` is a **GRCh38 1-based variant position**. It is not used for the
GRCh37 analysis and must not be compared directly with the GRCh37 BIM. Do not
subtract one from either position when forming a variant identifier. Liftover
assembly choice and coordinate-base conversion are separate concerns.

Audit identity: SHA-256
`7f45043f748cc6fad55fe4bc72e72151943896ce724ef2765f3add5eb763db6d`.

### LD reference

`LD_ref/g1k_EUR404_LDAK.{bed,bim,fam}` is the frozen European GRCh37 reference.
The BIM base-pair field is a **1-based variant position** and its variant IDs are
`chromosome:position`. The audited BIM contains 8,519,089 rows; SHA-256 is
`babd963bb32f1cc21b3c28e34456374218d7f93b40eac24291a254056a2713f2`.

## Current BEEA v6/v7 regulatory and annotation datasets

| Dataset | Assembly / coordinate contract | Role and coding decision |
| --- | --- | --- |
| JEME encoderoadmap lasso RDS files | GRCh37/hg19; `promoterFull` `location` is a 1-based TSS, while enhancer `START/END` parsed from strings such as `chrX:1008600-1009000` are already 0-based, half-open BED intervals | All tissues define the annotation/harmonization universe; E094 defines the Gastric network and HiC-overlap reference. Preserve exact enhancer node strings. Do not decrement JEME enhancer starts when writing LDAK gene files. The active v6/v7 LDAK promoter annotation does not use JEME `location`; it rematches the harmonized symbol to NCBI. |
| Jung HiC promoter-other RDS | GRCh37/hg19; promoter is a gene-symbol label with no promoter coordinate; `Interacting_fragment` is a 0-based, half-open BED-style HindIII-fragment interval | Gastric PO contacts contribute enhancer-promoter edges after selected-JEME overlap filtering. Active LDAK promoter coordinates are obtained by gene-symbol/alias matching to NCBI, not from the HiC workbook. Other-end fragment bounds pass unchanged to LDAK and bedtools. |
| Jung HiC promoter-promoter RDS | GRCh37; no genomic interval columns in the cached table | Gastric PP contributes promoter-promoter edges. PP must not be subjected to enhancer-coordinate filters. |
| NCBI GRCh37.p13 feature table | Source gene `start/end` are 1-based inclusive | Supplies promoter gene bounds and strand. In 0.1.1, `+` TSS is source `START`, `-` TSS is source `END`, a one-base promoter is `[TSS-1,TSS)`, and a gene body is `[START-1,END)`. |
| GENCODE v39 GRCh38 cleaned RDS | GRCh38; GTF-derived `start/end` are 1-based inclusive | Used for ENSG/gene-symbol harmonization only. Its GRCh38 coordinates must never be substituted for GRCh37 LDAK promoter coordinates. |
| Legacy `alias_link_nodup.rds` | No genomic coordinate contract | Maps symbols/aliases for frozen v6/v7 behavior. The unmarked legacy mapping direction must be preserved; do not silently replace it with the curated package default. |
| `ARACNe1600/GEjunc/network.txt` | No genomic coordinates; symbol/ENSG endpoints | Supplies GEjunc regulator-target edges with MI and p-value. Coordinate conversion is inapplicable. Preserve endpoint strings until explicit harmonization. |

The cached JEME and HiC upstream source formats do not themselves carry a
machine-readable coordinate-system attribute. Their fields must not be assigned
one blanket convention: JEME combines a 1-based TSS with BED enhancer
intervals, while HiC combines symbol-only promoters with BED-style
restriction-fragment intervals. A future resource builder should attach assembly and field-level
coordinate-system metadata and validate it on load rather than infer it from
column names.

### JEME promoter and enhancer audit

The cached E094 lasso object
`.cache/JEME/lasso/encoderoadmap_lasso.92.rds` contains 12,369 interactions and
6,965 unique `(promoter, chromosome, location, strand)` tuples. Its
`promoterFull` field has the original JEME form
`ENSG$symbol$chromosome$location$strand`; `location` is the TSS.

Matching those tuples by symbol, chromosome and strand to gene records in the
cached NCBI GRCh37.p13 feature table yielded 6,245 unique matched JEME tuples.
Across the 7,172 joined gene rows, 739 JEME locations exactly equal the
strand-aware NCBI boundary (`START` on `+`, `END` on `-`), compared with only 43
at boundary minus one and 43 at boundary plus one. Exact matches occur on both
strands (375 `+`, 364 `-`). This direct evidence agrees with Ensembl's documented
1-based inclusive convention: JEME `location` is a **1-based hg19 TSS**.

The original JEME source repository describes the enhancer inputs with names
such as `*.csv.bed.overlap.*` and gives examples such as
`chrX:1008600-1009000`. These are BED-derived enhancer intervals, so their
parsed bounds are **0-based, half-open**. This differs intentionally from the
TSS embedded in the same JEME record.

Audit identities:

| Resource | SHA-256 / revision |
| --- | --- |
| cached E094 JEME RDS | `511b03b1406136066ad56ef138b37e89f597dedaba3a784f45ead83751aa466b` |
| cached NCBI GRCh37.p13 feature table | `de711f7133885e4ada533ac2dea3c6baf5caaba77827a2138c598342dd892bc9` |
| original JEME source repository inspected | commit `f3c003a242d1e7b084d4a39397c36aed0e745c96` |

### HiC promoter-coordinate path

The cached Jung PO table contains `Promoter` gene symbols and separate
`Interacting_fragment` strings; the PP table contains only promoter symbols.
For the active v6/v7 LDAK annotation path, hotnetR2 harmonizes those symbols and
joins them to the NCBI GRCh37.p13 feature table. Consequently, the HiC promoter
gene bounds entering TSS construction are **NCBI 1-based inclusive** coordinates.
The same NCBI rematch is also used for JEME promoter symbols, so raw JEME TSS
locations do not control the active LDAK promoter `.loc` file.

This conclusion applies to promoter coordinates only. HiC PO
`Interacting_fragment` bounds describe HindIII restriction fragments and are
copied from the Jung supplementary workbook; they are not generated by matching
gene symbols to NCBI. Jung's fragment construction takes midpoints of adjacent
HindIII sites and merges nearby fragments. Local hg19 sequence spot checks of
the cached endpoints agree with those midpoint boundaries. hotnetR2 therefore
declares these fragment intervals BED-style and does not shift their starts.

## Coding decisions

1. Construct GWAS/LDAK predictor IDs from 1-based variant positions exactly as
   `chromosome:position`; never apply a BED `-1` conversion.
2. Match GWAS variants to the BIM by exact predictor first, then perform explicit
   direct/swapped/complement/complement-swapped allele orientation. Coordinate
   and allele harmonization must remain separate audit steps.
3. Treat source genomic intervals and LDAK/BED intervals as distinct types in
   code and documentation. Convert a 1-based inclusive `[s,e]` to `[s-1,e)`
   exactly once; do not convert an interval that is already BED, such as a JEME
   enhancer.
4. For canonical NCBI TSS conversion use source `START` for `+`, source `END`
   for `-`, then write `[TSS-1,TSS)`. Version 0.1.1 applies this correction;
   older v6/v7 result directories retain the superseded arithmetic.
5. Do not use GENCODE v39/GRCh38 positions in GRCh37 LDAK annotations. The file
   is an identifier mapping resource in these profiles.
6. Preserve exact enhancer and predictor strings in output and comparison code.
   Parsing a node label does not authorize normalizing or shifting its position.
7. Keep assembly, coordinate convention, source hash and conversion history in
   future immutable run/resource manifests.
8. Treat coordinate conventions at field level. In particular, JEME
   `location` is 1-based but JEME enhancer `START/END` are already BED; HiC
   promoter coordinates come from the NCBI rematch, whereas HiC other-end
   fragment coordinates are already BED.

## Scope and unresolved provenance

Raw BEEA processing and raw GEjunc ARACNe generation remain deferred upstream
TODOs. Other retained GWAS cohorts and historical analysis profiles have not
been coordinate-audited by this document and must not inherit the BEEA/MVP
conclusion solely because their columns have similar names. This reference must
be extended with dataset-specific evidence before those workflows are changed.
Corrected runs must use new output directories rather than overwriting the
historical v6/v7 artifacts generated with the superseded promoter/classifier
arithmetic.
