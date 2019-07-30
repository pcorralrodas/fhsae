//Moldova SAE - data preparation
//Sandu Cojocaru, Paul Corral
set more off
clear all

//global main "C:\Users\Utilizator\Desktop\10. MDA SAE\"
if (upper("`c(username)'")=="WB378870") global main   "C:\Users\WB378870\OneDrive - WBG\10. MDA SAE\"
else global main   "C:\Users\\`c(username)'\WBG\Paul Andres Corral Rodas - 10. MDA SAE\"
global hbs    "$main\1.data\HBS\"
global output "$main\1.data\output\"

*===============================================================================
// Open data, and include location means
*===============================================================================
use "$output\census territ means.dta", clear

	drop c_* r_*
	
	rename sh_* c_sh_*
	rename d_* c_d_*
	
	/*
	foreach x of varlist c_*{
		qui: sum `x'
		if ((r(sd)/r(mean))<0.4) drop `x'
	}	
	*/
tempfile loc_mu
save `loc_mu'

use "$output\hhlev_1cl_2p_10perc.dta", clear
preserve
gen hhid_n = _n
merge m:1 hbs using `loc_mu'
	drop if _m==2
	drop _m
	
	gen raion2 = int(raion/1e3)
	replace raion2 = int(raion2/1e1) if raion2>=1000
	
		
	gen hid2 = string(hbs)+"0" if length(string(hbs))==3
	replace hid2 = string(hbs) if hid2==""
	
	gen _HID = string(zona)+string(raion2)+hid2
	
	gen double HID = real(_HID)
	
	gen lnhh = ln(hh_size)
		
tempfile mycensus
save `mycensus'
restore

	groupfunction, first(c_* r_*) norestore by(hbs)

tempfile candr
save `candr'

//GET ID
use `mycensus'
	groupfunction, first(raion2 zona) by(HID hbs)
	
tempfile ids
save `ids'

use "$output\readyforfh.dta", clear

	merge 1:1 HID using `ids'
		drop _m
	
	merge 1:1 hbs using `loc_mu'
		drop _m
	drop raion
	gen raion = int(HID/1e4)	
	sort HID
	saveold "$output\fhR.dta", replace version(11)

		
local tot = 0

foreach x of varlist c_*{
	local myvars `myvars' `x'
}

/*
global myX c_sh_male c_d_student c_d_main_own c_d_self_ag_a
fhsae fgt0_hid $myX , re(fgt0_var_hid) fh(fh_fgt0) fhse(fh_fgt0_se) ///
fhcv(fh_fgt0_cv) gamma(fh_fgt0_gamma) out noneg method(reml)
cap rename hid HID
keep HID  fh_fgt0 fh_fgt0_se fh_fgt0_gamma fgt0_hid fgt0_var_hid population newpopw popw*

cap drop raion
gen raion = int(HID/1e4)

preserve
	groupfunction, sum(population) by(raion)
	putmata Naggr = population if raion!=.	
	putmata raion = raion if raion!=.
restore

levelsof raion, local(_mr)

local myraion
foreach x of local _mr{
	gen _`x' = raion==`x' if raion!=.
	local myraion `myraion' _`x'
}

putmata R = (`myraion')  if raion!=.
putmata narea = population if raion!=.
putmata mse = fh_fgt0_se if raion!=.
mata: mse = mse:^2
mata: RN=narea:*R
putmata Y = fh_fgt0 if raion!=.


mata: myaggr_mse = quadcross(RN,(narea:*mse)):/(Naggr:^2)
mata: myaggr     = quadcross(RN,Y):/Naggr
mata: D=(raion,myaggr,myaggr_mse)


clear 
getmata (HID fgt0 mse)=D

gen source = "EBLUP FH Stata"

append using "$output\raionFH_R.dta"
replace source = "FH R" if missing(source)

*/



global myX c_sh_male c_d_student c_d_main_own c_d_self_ag_a
cap run "C:\Users\WB378870\OneDrive - WBG\000.my_ados\fhsae\fhsae.ado"
cap drop raion
gen raion = int(HID/1e4)
/*
fhsae fgt0_hid $myX, re(fgt0_var_hid) fh(fh_fgt0) fhse(fh_fgt0_se) ///
fhcv(fh_fgt0_cv) gamma(fh_fgt0_gamma) out noneg method(reml)
sss
*/
//set trace on
//set traced 1
fhsae fgt0_hid $myX , re(fgt0_var_hid) fh(fh_fgt0) fhse(fh_fgt0_se) ///
fhcv(fh_fgt0_cv) gamma(fh_fgt0_gamma) out noneg method(reml) censuspop(population) ///
aggarea(raion) force
