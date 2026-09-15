#' JEME Tissue Mapping Key
#'
#' A dataset containing mapping information for JEME (Joint Effects of Mutations and Enhancers) results.
#' This key helps associate tissue numeric identifiers with descriptive names and categories.
#'
#' @format A data frame with information on JEME tissues:
#' \describe{
#'   \item{nfile}{Numeric file identifier used in JEME filenames}
#'   \item{file}{Filename of the JEME dataset}
#'   \item{cat1}{Main category/system}
#'   \item{cat2}{Sub-category/specific tissue}
#'   \item{desc1}{Long description of the tissue}
#'   \item{desc2}{Additional description or source info}
#'   \item{Epigenome mnemonic}{Roadmap epigenome mnemonic}
#'   \item{input}{Input identifier}
#'   \item{input_name}{Input name}
#'   \item{Sample group}{Sample grouping}
#'   \item{Sample name}{Sample description}
#'   \item{tiss}{Tissue identifier}
#' }
#' @source \url{https://labs.sbpdiscovery.org/centerandlabs/cancercenter/YipLab/jeme/}
"key_tissues_jeme"

#' Human Transcription Factor List (1600+)
#'
#' A comprehensive list of human transcription factors based on the Cell 2018 paper.
#'
#' @format A data frame with 1639 rows and 32 variables. `gene.ID` is the
#' Ensembl gene identifier, `gene.Name` the gene symbol, `gene.DBD` the
#' DNA-binding-domain class, and `Is.TF?` the transcription-factor assessment.
#' The remaining columns retain the source's curator assessments (`c1.`, `c2.`),
#' disagreements (`cu.`), final assessments (`f.`), curator notes (`n.`), and
#' external database annotations (`pc.`).
#' @source \url{https://doi.org/10.1016/j.cell.2018.01.029}
"tf_human_1600"

#' Selection of 25 Transcription Factors
#'
#' A prioritized list of 25 transcription factors found in `data-raw/TF/TF25.tsv`.
#'
#' @format A data frame with one column:
#' \describe{
#'   \item{gene_symbol}{HUGO Gene Symbol}
#' }
"tf_25"

#' Curated gene alias mappings
#'
#' September 2026 JEME and HiC workbook mappings override covered legacy
#' identifiers; uncovered legacy identifiers are retained. C11orf48 maps to
#' LBHD1 and MEGT1 to LY6G6D by explicit project curation. These are frozen
#' assets, not live HGNC lookups. The mapping_direction attribute is
#' raw_to_approved: gene is the original label and Hsym is the final symbol.
#' @format Data frames with 13 columns: merge status `_mrg`, source identifiers
#' `gsm`, `hh`, `ensg`, `gene`, HGNC-related fields `Hsym`, `Hhgnc_id`,
#' `hgnc_id`, `Hentrez_id`, `Hname`, `Hlocation`, and type fields `ty`, `Hty`.
#' `alias` contains 60422 rows; `alias_nodup` contains 59241 unique gene keys.
#' Unavailable metadata for newly added keys is NA. Source and conflict audits
#' and the external R builder are distributed with the installed package.
#' @name alias
#' @aliases alias_nodup alias_link alias_link_nodup
#' @usage data(alias_link)
#' data(alias_link_nodup)
#' @docType data
NULL

#' Frozen gene annotation table
#'
#' Gene annotation asset inherited from hotnetR for compatibility.
#' @format A data frame with 57853 rows and 20 columns retaining genomic
#' coordinates, strand, GTF feature fields, gene and transcript identifiers,
#' gene names and types, HGNC identifiers and annotation tags.
#' @name gencode_v39_genes_nodup
#' @aliases gtf
#' @usage data(gencode_v39_genes_nodup)
#' @docType data
NULL
