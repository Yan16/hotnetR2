*HGNC ~ direct download
/*
*https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
cd c:\_mb
clear
insheet using 2026_0907_hgnc_complete_set.txt, names
foreach var of varlist _all {
	rename `var' H`var'
}
rename Hsym symbol
codebook sym		//Unique values: 45,045                    Missing "": 0/45,045
codebook Hprev_sym  //Unique values: 12,658                    Missing "": 32,342/45,045  ~multiple ids per row(!)
codebook Hens	    //Unique values: 42,358                    Missing "": 2,684/45,045
gen priority=1
order priority
sort symbol
save 2026_0907_hgnc_complete_set_SORTsym.dta, replace
*/

*prev_sym
cd c:\_mb
clear
use 2026_0907_hgnc_complete_set_SORTsym.dta
keep Hhgnc_id symbol Hname Hlocus_group Hstat Hlocation Halias_sym Hprev_sym Hgene_group Hens
* 1. Split the variable Promoter into two new columns based on the semicolon
split Hprev_symbol, parse("|") gen(Hprev_symbol_)
* 2. Create a unique identifier for each original row so reshape works
gen id = _n
* 3. Reshape the data from wide to long format to create separate rows
reshape long Hprev_symbol_, i(id) j(sub_id)
* 4. Clean up by dropping empty rows (if any) and the temporary id variables
drop if missing(Hprev_symbol_)
drop id sub_id Hprev_symbol
* 5. Rename the new variable back to the original name
rename Hprev_symbol_ Hprev_symbol
gen id=Hprev_sym
replace Hprev_sym="id=prev"
replace Halias_sym="."
gen priority=2
order id priority
sort id
save 2026_0907_HGNC_prevRS_SORTprev_sym.dta, replace
codebook id		   //Unique values: 15,746                    Missing "": 0/15,882
codebook Hprev_sym //Unique values: 15,746                    Missing "": 0/15,882
codebook sym       //Unique values: 12,703                    Missing "": 0/15,882
codebook Hens	   //Unique values: 11,858                    Missing "": 1,056/15,882

*alias_sym
cd c:\_mb
clear
use 2026_0907_hgnc_complete_set_SORTsym.dta
keep Hhgnc_id symbol Hname Hlocus_group Hstat Hlocation Halias_sym Hprev_sym Hgene_group Hens
* 1. Split the variable Promoter into two new columns based on the semicolon
split Halias_symbol, parse("|") gen(Halias_symbol_)
* 2. Create a unique identifier for each original row so reshape works
gen id = _n
* 3. Reshape the data from wide to long format to create separate rows
reshape long Halias_symbol_, i(id) j(sub_id)
* 4. Clean up by dropping empty rows (if any) and the temporary id variables
drop if missing(Halias_symbol_)
drop id sub_id Halias_symbol
* 5. Rename the new variable back to the original name
rename Halias_symbol_ Halias_symbol
gen id=Halias_sym
replace Halias_sym="id=alias"
replace Hprev_sym="."
gen priority=3
order id priority
sort id
save 2026_0907_HGNC_aliasRS_SORTalias_sym.dta, replace
codebook id			//Unique values: 43,406                    Missing "": 0/44,787
codebook Halias_sym //1
codebook sym 	    //Unique values: 22,642                    Missing "": 0/44,787
codebook Hens		//Unique values: 21,576                    Missing "": 1,202/44,787

*MERGED (long)
clear
use 2026_0907_hgnc_complete_set_SORTsym.dta
keep priority Hhgnc_id symbol Hname Hlocus_group Hstat Hlocation Halias_sym Hprev_sym Hgene_group Hens
gen id=symbol
order id
codebook id			//Unique values: 45,045                    Missing "": 0/45,045
codebook sym		//Unique values: 45,045                    Missing "": 0/45,045
codebook Hprev_sym  //Unique values: 12,658                    Missing "": 32,342/45,045 *multi/row
codebook Halias_sym //Unique values: 22,572                    Missing "": 22,403/45,045 *multi/row
codebook Hens	    //Unique values: 42,358                    Missing "": 2,684/45,045
sort id
append using 2026_0907_HGNC_prevRS_SORTprev_sym.dta
append using 2026_0907_HGNC_aliasRS_SORTalias_sym.dta
save 2026_0907_HGNC_APPENDlong_SORTgene.dta, replace
export excel using 2026_0907_HGNC_APPENDlong_SORTgene.xlsx, firstrow(var) replace
codebook id			//Unique values: 103,137                   Missing "": 0/105,714
codebook sym		//Unique values: 45,045                    Missing "": 0/105,714
codebook Hprev_sym  //Unique values: 12,660                    Missing "": 32,342/105,714
codebook Halias_s	//Unique values: 22,574                    Missing "": 22,403/105,714
codebook Hens	    //Unique values: 42,358                    Missing "": 4,942/105,714


/*
The HGNC now archive the complete HGNC dataset file (both tab separated and JSON formats) for each month and each quarter and store them in our Google Cloud Storage Bucket.
File fields/columns information

The monthly files are produced at the same time as the first update of every month while the quarterly files are produced at the same time as the first update of Jan, Apr, Jul & Oct. If a monthly file is over 365 days old, the file will be deleted for us to save disk space however, the quarterly files are not deleted currently so for snapshots older than 365 days please use the quarterly files.

There are essentially two types of data file (excluding the file format type) of hgnc_complete_set and withdrawn. The hgnc_complete_set is a set of all approved gene symbol reports found on the GRCh38 reference and the alternative reference loci (see fig. 1 for a list of columns/headings).

hgnc_id                  = HGNC ID. A unique ID created by the HGNC for every
                           approved symbol. 

symbol                   = The HGNC approved gene symbol. Equates to the
                           "APPROVED SYMBOL" field within the gene symbol
                           report.

name                     = HGNC approved name for the gene. Equates to the
                           "APPROVED NAME" field within the gene symbol report.

locus_group              = A group name for a set of related locus types as
                           defined by the HGNC (e.g. non-coding RNA).

locus_type               = The locus type as defined by the HGNC (e.g. RNA,
                           transfer).

status                   = Status of the symbol report, which can be either
                           "Approved" or "Entry Withdrawn".

location                 = Cytogenetic location of the gene (e.g. 2q34).

location_sortable        = Same as "location" but single digit chromosomes are
                           prefixed with a 0 enabling them to be sorted in
                           correct numerical order (e.g. 02q34).

alias_symbol             = Other symbols used to refer to this gene as seen in
                           the "SYNONYMS" field in the symbol report. 

alias_name               = Other names used to refer to this gene as seen in
                           the "SYNONYMS" field in the gene symbol report.

prev_symbol              = Symbols previously approved by the HGNC for this
                           gene. Equates to the "PREVIOUS SYMBOLS & NAMES" field
                           within the gene symbol report.

prev_name                = Gene names previously approved by the HGNC for this
                           gene. Equates to the "PREVIOUS SYMBOLS & NAMES" field
                           within the gene symbol report.

gene_family              = Name given to a gene family or group the gene has been
                           assigned to. Equates to the "GENE FAMILY" field within
                           the gene symbol report.

gene_family_id           = ID used to designate a gene family or group the gene
                           has been assigned to.

date_approved_reserved   = The date the entry was first approved.

date_symbol_changed      = The date the gene symbol was last changed.

date_name_changed        = The date the gene name was last changed.

date_modified            = Date the entry was last modified.

entrez_id                = Entrez gene ID. Found within the "GENE RESOURCES"
                           section of the gene symbol report.

ensembl_gene_id          = Ensembl gene ID. Found within the "GENE RESOURCES"
                           section of the gene symbol report.

vega_id                  = Vega gene ID. Found within the "GENE RESOURCES"
                           section of the gene symbol report.

ucsc_id                  = UCSC gene ID. Found within the "GENE RESOURCES"
                           section of the gene symbol report.

ena                      = International Nucleotide Sequence Database
                           Collaboration (GenBank, ENA and DDBJ) accession
                           number(s). Found within the "NUCLEOTIDE SEQUENCES"
                           section of the gene symbol report.

refseq_accession         = RefSeq nucleotide accession(s). Found within the
                           "NUCLEOTIDE SEQUENCES" section of the gene symbol
                           report.

ccds_id                  = Consensus CDS ID. Found within the
                           "NUCLEOTIDE SEQUENCES" section of the gene symbol
                           report.

uniprot_ids              = UniProt protein accession. Found within the
                           "PROTEIN RESOURCES" section of the gene symbol
                           report.

pubmed_id                = Pubmed and Europe Pubmed Central PMID(s).

mgd_id                   = Mouse genome informatics database ID. Found within
                           the "HOMOLOGS" section of the gene symbol report.

rgd_id                   = Rat genome database gene ID. Found within the
                           "HOMOLOGS" section of the gene symbol report.

lsdb                     = The name of the Locus Specific Mutation Database and
                           URL for the gene separated by a | character

cosmic                   = Symbol used within the Catalogue of somatic
                           mutations in cancer for the gene.

omim_id                  = Online Mendelian Inheritance in Man (OMIM) ID

mirbase                  = miRBase ID

homeodb                  = Homeobox Database ID

snornabase               = snoRNABase ID

bioparadigms_slc         = Symbol used to link to the SLC tables database at
                           bioparadigms.org for the gene

orphanet                 = Orphanet ID

pseudogene.org           = Pseudogene.org

horde_id                 = Symbol used within HORDE for the gene

merops                   = ID used to link to the MEROPS peptidase database

imgt                     = Symbol used within international ImMunoGeneTics
                           information system

iuphar                   = The objectId used to link to the IUPHAR/BPS Guide to
                           PHARMACOLOGY database. To link to IUPHAR/BPS Guide
                           to PHARMACOLOGY database only use the number
                           (only use 1 from the result objectId:1)

kznf_gene_catalog        = ID used to link to the Human KZNF Gene Catalog

mamit-trnadb             = ID to link to the Mamit-tRNA database

cd                       = Symbol used within the Human Cell Differentiation
                           Molecule database for the gene

lncrnadb                 = lncRNA Database ID - Resource now defunct

enzyme_id                = ENZYME EC accession number

intermediate_filament_db = ID used to link to the Human Intermediate Filament
                           Database

agr                      = The HGNC ID that the Alliance of Genome Resources
                           (AGR) have linked to their record of the gene. Use
                           the HGNC ID to link to the AGR.

mane_select              = NCBI and Ensembl transcript IDs/acessions
                           including the version number for one high-quality
                           representative transcript per protein-coding gene
                           that is well-supported by experimental data and
                           represents the biology of the gene. The IDs are
                           delimited by |.

Figure 1: Columns/headings within the hgnc_complete_set files.

The withdrawn file contains all gene symbol reports that are no longer approved. Either the symbol has been withdrawn or merged/split into another report (see fig. 2 for a list of columns/headings).

HGNC_ID               = The HGNC ID of the withdrawn record.
STATUS                = Can either be "Entry Withdrawn" or "Merged/Split"
WITHDRAWN_SYMBOL      = The symbol of the withdrawn record.
MERGED_INTO_REPORT(S) = Shows what record(s) replaced the withdrawn record if the
                        status is "Merged/Split". Each replacement has the format
                        HGNC_ID|SYMBOL|STATUS. If the withdrawn record is split,
                        their will be more than one replacement and they will be
                        comma separated.

Figure 2: Columns/headings within the withdrawn files.
*/
