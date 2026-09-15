*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2026-0905
*RAW input [HiC]
/*
*HiC P-O
*https://www.nature.com/articles/s41588-019-0494-8#Sec23
cd "\\vpensbst\bstshared\Epidemiology\Buas\GRANTS\R01_fxGWAS\!draft_chr04_LINC\LY_hotnet\RAW_jeme_hiC"
clear
import excel using 41588_2019_494_MOESM3_ESM_m1.xlsx, firstrow  //enh
sort Pro
save hic_po_input.dta, replace

*HiC P-P
cd "\\vpensbst\bstshared\Epidemiology\Buas\GRANTS\R01_fxGWAS\!draft_chr04_LINC\LY_hotnet\RAW_jeme_hiC"
clear
import excel using 41588_2019_494_MOESM4_ESM_m1.xlsx, firstrow //pro
save hic_pp_input.dta, replace
*/

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ hic.po> split multi-gene entries; update excel gene_name typos; 
cd "\\vpensbst\bstshared\Epidemiology\Buas\GRANTS\R01_fxGWAS\!draft_chr04_LINC\LY_hotnet\RAW_jeme_hiC"
clear
use hic_po_input.dta

* 1. Split the variable Promoter into two new columns based on the semicolon
split Promoter, parse(";") gen(Promoter_)
* 2. Create a unique identifier for each original row so reshape works
gen id = _n
* 3. Reshape the data from wide to long format to create separate rows
reshape long Promoter_, i(id) j(sub_id)
* 4. Clean up by dropping empty rows (if any) and the temporary id variables
drop if missing(Promoter_)
drop id sub_id Promoter
* 5. Rename the new variable back to the original name
rename Promoter_ Promoter
gen c=substr(Int,1,strpos(Int,".")-1)
sort Pro
ta c Pro in 1/856

/*
DEC1	DELEC1	9q33.1	Https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:23658		YES
DEC1	BHLHE40 3p26.1	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1046		X

MAR1	RTL1 	14q32.2	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:14665		X
MARC1	MTARC1	1q41 	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:26189		YES (N=9)
MARCH1	MARCHF1	4q32.2	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:26077		YES (N=117)

SEP1	XRN1	3q23 	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:30654		X
SEPT1	SEPTIN1	16p11.2	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:2879		YES

MAR2	PEG10	7q21.3	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:14005
MARC2	MTARC2	1q41	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:26064		YES(N=29)
MARCH2	MARCHF2	19p13.2 https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:28038		YES(N=5)

MAR3	Xq21.1 
MARCH3	MARCHF3	5q23.2	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:28728		YES

SEP3	SEPTIN3	22q13.2	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:10750		YES
SEPT3	~same

MAR4	Xq23
MARCH4	MARCHF4	2q35	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:29269		YES

MAR6	22q13.31 
MARCH6	MARCHF6	5p15.2 	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:30550		YES

SEPT6	SEPTIN6	Xq24 	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:15848		YES

MAR7	Xq27.1 
MARCH7	MARCHF7	2q24.2 	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:17393		YES

SEPT7	SEPTIN7	7p14.2	https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1717		YES

MAR8	Xq26.3
MARCH8	MARCHF8	10q11.21 https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:23356	YES

SEPT9	SEPTIN9	17q25.3	 https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:7323		YES

MARCH10	MARCHF10 17q23.2 https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:26655	YES

MARCH11	MARCHF11 5p15.1	 https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:33609	YES

SEPT11	SEPTIN11 4q21.1  https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:25589	YES

           |                                                              Promoter
         c |     1-Dec      1-Mar XX   1-Sep      2-Mar XX   3-Mar      3-Sep      4-Mar      6-Mar      6-Sep      7-Mar      7-Sep      8-Mar |     Total
-----------+------------------------------------------------------------------------------------------------------------------------------------+----------
      chr1 |         0          9 X        0         29          0          0          0          0          0          0          0          0 |        38 
     chr10 |         0          0          0          0          0          0          0          0          0          0          0         28 |        28 
     chr16 |         0          0        152          0          0          0          0          0          0          0          0          0 |       152 
     chr17 |         0          0          0          0          0          0          0          0          0          0          0          0 |       170 
     chr19 |         0          0          0          5  X       0          0          0          0          0          0          0          0 |         5 
      chr2 |         0          0          0          0          0          0          4          0          0          7          0          0 |        11 
     chr22 |         0          0          0          0          0        164          0          0          0          0          0          0 |       164 
      chr4 |         0        117          0          0          0          0          0          0          0          0          0          0 |       151 
      chr5 |         0          0          0          0         20          0          0         47          0          0          0          0 |        73 
      chr7 |         0          0          0          0          0          0          0          0          0          0         30          0 |        30 
      chr9 |        30          0          0          0          0          0          0          0          0          0          0          0 |        30 
      chrX |         0          0          0          0          0          0          0          0          4          0          0          0 |         4 
-----------+------------------------------------------------------------------------------------------------------------------------------------+----------
     Total |        30        126        152         34         20        164          4         47          4          7         30         28 |       856 

           |                  Promoter
         c |     9-Sep     10-Mar     11-Mar     11-Sep |     Total
-----------+--------------------------------------------+----------
      chr1 |         0          0          0          0 |        38 
     chr10 |         0          0          0          0 |        28 
     chr16 |         0          0          0          0 |       152 
     chr17 |       155         15          0          0 |       170 
     chr19 |         0          0          0          0 |         5 
      chr2 |         0          0          0          0 |        11 
     chr22 |         0          0          0          0 |       164 
      chr4 |         0          0          0         34 |       151 
      chr5 |         0          0          6          0 |        73 
      chr7 |         0          0          0          0 |        30 
      chr9 |         0          0          0          0 |        30 
      chrX |         0          0          0          0 |         4 
-----------+--------------------------------------------+----------
     Total |       155         15          6         34 |       856 
*/

codebook Pro // Unique values: 18,652                    Missing "": 0/932,357
replace Promoter=subinstr(Promoter," ","",.)
replace Promoter="DELEC1" if Promoter=="1-Dec" & c=="chr9"
replace Promoter="MTARC1" if Promoter=="1-Mar" & c=="chr1"
replace Promoter="MARCHF1" if Promoter=="1-Mar" & c=="chr4"
replace Promoter="SEPTIN1" if Promoter=="1-Sep" & c=="chr16"
replace Promoter="MTARC2" if Promoter=="2-Mar" & c=="chr1"
replace Promoter="MARCHF2" if Promoter=="2-Mar" & c=="chr19"
replace Promoter="MARCHF3" if Promoter=="3-Mar" & c=="chr5"
replace Promoter="SEPTIN3" if Promoter=="3-Sep" & c=="chr22"
replace Promoter="MARCHF4" if Promoter=="4-Mar" & c=="chr2"
replace Promoter="MARCHF6" if Promoter=="6-Mar" & c=="chr5"
replace Promoter="SEPTIN6" if Promoter=="6-Sep" & c=="chrX"
replace Promoter="MARCHF7" if Promoter=="7-Mar" & c=="chr2"
replace Promoter="SEPTIN7" if Promoter=="7-Sep" & c=="chr7"
replace Promoter="MARCHF8" if Promoter=="8-Mar" & c=="chr10"
replace Promoter="SEPTIN9" if Promoter=="9-Sep" & c=="chr17"
replace Promoter="MARCHF10" if Promoter=="10-Mar" & c=="chr17"
replace Promoter="MARCHF11" if Promoter=="11-Mar" & c=="chr5"
replace Promoter="SEPTIN11" if Promoter=="11-Sep" & c=="chr4"
codebook Pro // Unique values: 18,654                    Missing "": 0/932,357
codebook Int // Unique values: 200,595                   Missing "": 0/932,357
ta c Pro in 1/856
drop c
save hic_po_input_uniqueF.dta, replace


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ hic.pp> split multi-gene entries; update excel gene_name typos; 
cd "\\vpensbst\bstshared\Epidemiology\Buas\GRANTS\R01_fxGWAS\!draft_chr04_LINC\LY_hotnet\RAW_jeme_hiC"
clear
use hic_pp_input.dta, replace
codebook Pro //Unique values: 9,853                     Missing "": 0/79,989
codebook B   //Unique values: 9,900                     Missing "": 0/79,989

replace Promoter=subinstr(Promoter," ","",.)
replace B=subinstr(B," ","",.)

codebook Pro //Unique values: 9,853                     Missing "": 0/79,989
codebook B   //Unique values: 9,900                     Missing "": 0/79,989

sort Pro
sort B
/*
*verify CHR for [paired gene] when multiple gene possibilities for Excel_typo_labels (2-Mar; 1-Mar)
Promoter	B	Tissue_type
 2-Mar	HNRNPM	Spleen					19p13.2
 2-Mar	MYO1F	Spleen					19p13.2 
 2-Mar	MYO1F	Neural Progenitor Cell
 2-Mar	HNRNPM	Psoas
 
 Promoter	B	Tissue_type
ANGPTL4	 2-Mar	Adrenal Gland			19p13.2 

"2-Mar" --> "MARCHF2" ~~~~for all cases above
*/
*
replace Promoter="DELEC1" if Promoter=="1-Dec"
*replace Promoter="MTARC1" if Promoter=="1-Mar"
replace Promoter="MARCHF1" if Promoter=="1-Mar"
replace Promoter="SEPTIN1" if Promoter=="1-Sep"
*replace Promoter="MTARC2" if Promoter=="2-Mar"
replace Promoter="MARCHF2" if Promoter=="2-Mar"
replace Promoter="MARCHF3" if Promoter=="3-Mar"
replace Promoter="SEPTIN3" if Promoter=="3-Sep"
replace Promoter="MARCHF4" if Promoter=="4-Mar"
replace Promoter="MARCHF6" if Promoter=="6-Mar"
replace Promoter="SEPTIN6" if Promoter=="6-Sep"
replace Promoter="MARCHF7" if Promoter=="7-Mar"
replace Promoter="SEPTIN7" if Promoter=="7-Sep"
replace Promoter="MARCHF8" if Promoter=="8-Mar"
replace Promoter="SEPTIN9" if Promoter=="9-Sep"
replace Promoter="MARCHF10" if Promoter=="10-Mar"
replace Promoter="MARCHF11" if Promoter=="11-Mar"
replace Promoter="SEPTIN11" if Promoter=="11-Sep"

replace B="DELEC1" if B=="1-Dec"
*replace B="MTARC1" if B=="1-Mar"
replace B="MARCHF1" if B=="1-Mar"
replace B="SEPTIN1" if B=="1-Sep"
*replace B="MTARC2" if B=="2-Mar"
replace B="MARCHF2" if B=="2-Mar"
replace B="MARCHF3" if B=="3-Mar"
replace B="SEPTIN3" if B=="3-Sep"
replace B="MARCHF4" if B=="4-Mar"
replace B="MARCHF6" if B=="6-Mar"
replace B="SEPTIN6" if B=="6-Sep"
replace B="MARCHF7" if B=="7-Mar"
replace B="SEPTIN7" if B=="7-Sep"
replace B="MARCHF8" if B=="8-Mar"
replace B="SEPTIN9" if B=="9-Sep"
replace B="MARCHF10" if B=="10-Mar"
replace B="MARCHF11" if B=="11-Mar"
replace B="SEPTIN11" if B=="11-Sep"

codebook Pro //Unique values: 9,853                     Missing "": 0/79,989
codebook B   //Unique values: 9,900                     Missing "": 0/79,989

* 1. Split the variable Promoter into two new columns based on the semicolon
split Promoter, parse(";") gen(Promoter_)
* 2. Create a unique identifier for each original row so reshape works
gen id = _n
* 3. Reshape the data from wide to long format to create separate rows
reshape long Promoter_, i(id) j(sub_id)
* 4. Clean up by dropping empty rows (if any) and the temporary id variables
drop if missing(Promoter_)
drop id sub_id Promoter
* 5. Rename the new variable back to the original name
rename Promoter_ Promoter

* 1. Split the variable Promoter into two new columns based on the semicolon
split B, parse(";") gen(B_)
* 2. Create a unique identifier for each original row so reshape works
gen id = _n
* 3. Reshape the data from wide to long format to create separate rows
reshape long B_, i(id) j(sub_id)
* 4. Clean up by dropping empty rows (if any) and the temporary id variables
drop if missing(B_)
drop id sub_id B
* 5. Rename the new variable back to the original name
rename B_ B
sort Pro
save hic_pp_input_uniqueF.dta, replace


**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ hiC~FINAL PO
*c:\_mb\HGNC_SORTsym
cd "\\vpensbst\bstshared\Epidemiology\Buas\GRANTS\R01_fxGWAS\!draft_chr04_LINC\LY_hotnet\RAW_jeme_hiC"
clear
use hic_po_input_uniqueF
codebook Int //Unique values: 200,595                   Missing "": 0/932,357
codebook Pro //Unique values: 18,654                    Missing "": 0/932,357
egen r=rank(_n),by(Pro)
keep if r==1
drop r
codebook Pro //Unique values: 18,654                    Missing "": 0/18,654
keep Pro Int
rename Pro id

sort id
joinby id using c:\_mb\2026_0907_HGNC_APPENDlong_SORTid, unmatched(master)
ta _merge
codebook id if _merge==1 //Unique values: 874                       Missing "": 0/874
codebook id if _merge==3 //Unique values: 17,780                    Missing "": 0/18,347
/*

                       _merge |      Freq.     Percent        Cum.
------------------------------+-----------------------------------
          only in master data |        874        4.55        4.55
both in master and using data |     18,347       95.45      100.00
------------------------------+-----------------------------------
                        Total |     19,221      100.00
*/
****
drop if _merge==1
drop _merge
****

ta prior,mi
gen ch=substr(Hlocation,1,strpos(Hlocation,"p")-1)
replace ch=substr(Hlocation,1,strpos(Hlocation,"q")-1) if ch==""
*
gen t=substr(Int,4,.)
gen ichr=substr(t,1,strpos(t,".")-1)
drop t

gen k=0
replace k=1 if id==symbol & priority==1 & ichr==ch //n=1 xx NAA38
replace k=2 if id!=symbol & priority==2 & ichr==ch
replace k=3 if id!=symbol & priority==3 & ichr==ch

egen a=rank(_n),by(id)
egen xa=max(a),by(id)
gsort -xa id a
order id xa a k prior ichr ch Hlocation sym Hprev Halias

gen k1 byte=(k==1)
egen xk1=max(k1),by(id)
gen k2 byte=(k==2)
egen xk2=max(k2),by(id)
gen k3 byte=(k==3)
egen xk3=max(k3),by(id)
order xk1 xk2 xk3

*IF MULTIPLE ROWS per id--> if there is k=1 entry for given id, drop other rows; etc k=2,k=3
****
drop if xa>1 & xk1==1 & k!=1
****
drop a xa
egen a=rank(_n),by(id)
egen xa=max(a),by(id)
gsort -xa id a
order id xa a k xk1 xk2 xk3 prior ichr ch Hlocation sym Hprev Halias

****
drop if xa>1 & xk2==1 & k!=2
****
drop a xa
egen a=rank(_n),by(id)
egen xa=max(a),by(id)
gsort -xa id a
order id xa a k xk1 xk2 xk3 prior ichr ch Hlocation sym Hprev Halias

****
drop if xa>1 & xk3==1 & k!=3
****
drop a xa
egen a=rank(_n),by(id)
egen xa=max(a),by(id)
gsort -xa id a
order id xa a k xk1 xk2 xk3 prior ichr ch Hlocation sym Hprev Halias

codebook id		// Unique values: 17,780                    Missing "": 0/17,782
codebook sym	// Unique values: 17,753                    Missing "": 0/17,782
drop xa a k xk1 xk2 xk3 k1 k2 k3 

*check dups
egen a=rank(_n),by(sym)
egen xa=max(a),by(sym)
gsort -xa sym a
order xa a sym
drop xa a

*ID=hiC(Pro) ---> symbol[=final]
order id
sort id
save final_harmonize_hiC_PO.dta, replace
keep id symbol priority Hlocation Hhgnc Hname  Hlocus_group Hstat Hgene_group Hens 
sort id
export excel using final_harmonize_hiC_PO.xlsx, firstrow(var) replace



**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ hiC~FINAL PP
clear
use hic_pp_input_uniqueF.dta
codebook Pro //Unique values: 11,805                    Missing "": 0/133,361
codebook B   //Unique values: 11,858                    Missing "": 0/133,361

rename Pro id
sort id
joinby id using final_harmonize_hiC_PO.dta, unmatched(master)
ta _merge
/*
                       _merge |      Freq.     Percent        Cum.
------------------------------+-----------------------------------
          only in master data |      6,101        4.57        4.57
both in master and using data |    127,282       95.43      100.00
------------------------------+-----------------------------------
                        Total |    133,383      100.00
*/
codebook id if _merge==1 //Unique values: 721                       Missing "": 0/6,101
codebook id if _merge==3 //Unique values: 11,084                    Missing "": 0/127,282
****
drop if _merge==1
drop _merge
****
rename id Pro
rename sym Pro_sym
keep Tiss Pro Pro_sym B

*
rename B id
sort id
joinby id using final_harmonize_hiC_PO.dta, unmatched(master)
ta _merge
/*
                       _merge |      Freq.     Percent        Cum.
------------------------------+-----------------------------------
          only in master data |      6,045        4.75        4.75
both in master and using data |    121,248       95.25      100.00
------------------------------+-----------------------------------
                        Total |    127,293      100.00
*/
codebook id if _merge==1 //Unique values: 679                       Missing "": 0/6,045
codebook id if _merge==3 //Unique values: 11,042                    Missing "": 0/121,248
****
drop if _merge==1
drop _merge
****
rename id B
rename sym B_sym
keep Tiss Pro B Pro_sym B_sym

sort Pro
save final_harmonize_hiC_PP_pairs.dta, replace
export excel using final_harmonize_hiC_PP_pairs.xlsx, firstrow(var) replace

