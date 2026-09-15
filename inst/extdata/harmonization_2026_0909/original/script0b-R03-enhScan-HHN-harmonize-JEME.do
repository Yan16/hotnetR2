*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2026-0906 JEME
*JEME
/*
cd "\\vpensbst\bstshared\Epidemiology\Buas\GRANTS\R01_fxGWAS\!draft_chr04_LINC\LY_hotnet\RAW_jeme_hiC"
clear
insheet using L_EN_jeme.bed, double nonames
*(10 vars, 6,133,653 obs)
rename v5 gene
gen Jensg=substr(v4,1,strpos(v4,".")-1)
rename v4 J1ensg
sort gene
save L_EN_jeme_merged.sortGENE.dta, replace
*/

*MERGE to new_master_HGNC
cd "\\vpensbst\bstshared\Epidemiology\Buas\GRANTS\R01_fxGWAS\!draft_chr04_LINC\LY_hotnet\RAW_jeme_hiC"
clear
use L_EN_jeme_merged.sortGENE.dta
codebook J1ensg		//Unique values: 19,195                    Missing "": 0/6,133,653
codebook gene		//Unique values: 19,172                    Missing "": 0/6,133,653
egen r=rank(_n),by(gene)
keep if r==1
drop r
codebook J1ensg		//Unique values: 19,172                    Missing "": 0/19,172
codebook gene		//Unique values: 19,172                    Missing "": 0/19,172
keep J1ensg Jensg gene v1
gen id=gene
order id
****
sort id
joinby id using c:\_mb\2026_0907_HGNC_APPENDlong_SORTid, unmatched(master)
****
ta _merge
codebook id if _merge==1 //Unique values: 480                       Missing "": 0/480
codebook id if _merge==3 //Unique values: 18,692                    Missing "": 0/19,270
/*
                       _merge |      Freq.     Percent        Cum.
------------------------------+-----------------------------------
          only in master data |        480        2.43        2.43
both in master and using data |     19,270       97.57      100.00
------------------------------+-----------------------------------
                        Total |     19,750      100.00
*/
****
drop if _merge==1
drop _merge
****

ta prior,mi
gen ch=substr(Hlocation,1,strpos(Hlocation,"p")-1)
replace ch=substr(Hlocation,1,strpos(Hlocation,"q")-1) if ch==""

gen k=0
replace k=1 if id==symbol & priority==1 & substr(v1,4,.)==ch //n=1 xx NAA38
replace k=2 if id!=symbol & priority==2 & substr(v1,4,.)==ch
replace k=3 if id!=symbol & priority==3 & substr(v1,4,.)==ch

egen a=rank(_n),by(id)
egen xa=max(a),by(id)
gsort -xa id a
order id xa a k prior v1 ch Hlocation sym Hprev Halias

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
order id xa a k xk1 xk2 xk3 prior v1 ch Hlocation sym Hprev Halias

****
drop if xa>1 & xk2==1 & k!=2
****
drop a xa
egen a=rank(_n),by(id)
egen xa=max(a),by(id)
gsort -xa id a
order id xa a k xk1 xk2 xk3 prior v1 ch Hlocation sym Hprev Halias

****
drop if xa>1 & xk3==1 & k!=3
****
drop a xa
egen a=rank(_n),by(id)
egen xa=max(a),by(id)
gsort -xa id a
order id xa a k xk1 xk2 xk3 prior v1 ch Hlocation sym Hprev Halias

codebook id		// Unique values: 18,692                    Missing "": 0/18,694
codebook sym	// Unique values: 18,645                    Missing "": 0/18,694
drop xa a k xk1 xk2 xk3 k1 k2 k3 

*check dups
egen a=rank(_n),by(sym)
egen xa=max(a),by(sym)
gsort -xa sym a
order xa a sym
drop xa a

*ID=jeme ---> symbol[=final]
order id
sort id
save final_harmonize_jemeF.dta, replace
keep id J1ensg symbol priority Hlocation Hhgnc Hname  Hlocus_group Hstat Hgene_group Hens
export excel using final_harmonize_jemeF.xlsx, firstrow(var) replace

