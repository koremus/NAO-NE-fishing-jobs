/******************************  METADATA *********************************
generateFigure_1_2.do
Script Description: This STATA script inputs and merges data for regression. 
Code developed for STATA v15.

Author:  Kimberly L. Oremus
Date: April 21, 2019

Data from: 
NMFS Annual Commercial Landings data
regional, yearly landings
https://www.st.nmfs.noaa.gov/commercial-fisheries/commercial-landings/annual-landings/index
see SI appendix, Table S1 for a list of species that were aggrgated by region and year
data located in input folder "NE-NMFS-total-71-17.csv" and "SA-NMFS-total-71-17.csv"
accessed on April 19, 2019

NCAR ClimateDataGuide
DJFM annual NAO index
https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based
data located in input folder "NAO_slp_djfm_1864_2017.csv"
accessed on January 10, 2019
**************************************************************************/

/******************************  HEADER **********************************/
clear all
set trace off
set more off
capture log close

//set local directory 
cd ..

global inputDir "input"
global resultsDir "tables"
global figuresDir "figures"
global tempDir "temp"

/******************* Fig. 1A: NAO ************************/
insheet using $inputDir/NAO_slp_djfm_1864_2017.csv, comma names clear
rename nao nao_slp_djfm		
label var nao_slp_djfm "Winter (DJFM) SLP-based NAO index"	
drop if year<1965
gen nao_pos=nao_slp_djfm
replace nao_pos=0 if nao_pos<0
twoway (area nao_slp_djfm year, fcolor(blue)) ///
(area nao_pos year, lcolor(white) fcolor(red) xlab(1965(10)2017) ///
legend(label(1 "negative phase") label(2 "positive phase")) ytitle("NAO index (std dev)") ///
ylab(,nogrid) plotregion(fcolor(white)) graphregion(fcolor(white)))
graph export $figuresDir/fig1a-NAO-65-17.pdf, as(pdf) replace
 
/********* Fig. 2: NAO effects on catch and revenue **********/

//LOAD AND MERGE
//regional, aggregate landings data with NAO

	//Winter NAO index (station-based sea-level pressure (SLP))
	insheet using $inputDir/nao_slp_djfm_1864_2017.csv, comma names clear
	rename nao nao_slp_djfm		
	label var nao_slp_djfm "Winter (DJFM) SLP-based NAO index"	
	save $tempDir/nao_slp, replace
	//New England collapsed landings data
	insheet using $inputDir/NE-NMFS-total-71-17.csv, comma names clear
	merge m:1 year using $tempDir/nao_slp
	drop _merge
	sort year
	drop if year < 1965
	save $tempDir/select-ne, replace
	//South Alantic collapsed landings data
	insheet using $inputDir/SA-NMFS-total-71-17.csv, comma names clear
	merge m:1 year using $tempDir/nao_slp
	drop _merge
	sort year
	drop if year < 1965
	save $tempDir/select-sa, replace


//NEW ENGLAND PLOTS

	use $tempDir/select-ne, replace

	gen time1=year-1871
	gen time2 = time1^2
	gen time3=time1^3
	gen time4=time1^4
	gen time5=time1^5

//FIG. 1A: NE CATCH PLOT
	gen l_catch=ln(main_mtons)
	drop if year<1965
	
	local outcome "total catch"
	local acro "NE"
	local region "New England"
	local order "5"
	tsset year	

	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order', robust bw(auto)
	gen sample=1 if e(sample)
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm
	est store L0
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm
	est store L1
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm
	est store L2
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm
	est store L3
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm
	est store L4
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm
	est store L5
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm+L6.nao_slp_djfm
	est store L6
	
	coefplot ///
	(L0\L1\L2\L3\L4\L5\L6, recast(line) lcolor(dkorange) lwidth(1) ciopts(recast(rarea) lcolor(dkorange%0) fcolor(dkorange%30) lwidth(1)) ///
	yline(0, lwidth(1) lc(gray) lp(dash)) ylabel(.15(.05)-.2, nogrid angle(0))), ///
	vertical aseq swapnames ///
	coeflabels(L0 = "0" L1="1" L2="2" L3="3" L4="4" L5="5" L6="6") ///
	ylab(,nogrid) graphregion(color(white)) ///
	ytitle("Cumulative effect of 1-unit NAO increase" "on percentage of regional `outcome'") ///
	xtitle(Annual lags) graphregion(color(white)) legend(label(2 "`region'")) 
	graph export $figuresDir/fig2a-catch-`acro'.pdf, as(pdf) replace

	drop sample _est_L0 _est_L1 _est_L2 _est_L3 _est_L4 _est_L5 _est_L6 l_catch

//FIG. 2B: NE REVENUE PLOT
	gen l_catch=ln(main_rev)
	local outcome "total revenue"

	tsset year

	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order', robust bw(auto)
	gen sample=1 if e(sample)
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm
	est store L0
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm
	est store L1
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm
	est store L2
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm
	est store L3
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm
	est store L4
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm
	est store L5
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm+L6.nao_slp_djfm
	est store L6

	coefplot ///
	(L0\L1\L2\L3\L4\L5\L6, recast(line) lcolor(dkorange) lwidth(1) ciopts(recast(rarea) lcolor(dkorange%0) fcolor(dkorange%30) lwidth(1)) ///
	yline(0, lwidth(1) lc(gray) lp(dash)) ylabel(.15(.05)-.2, nogrid angle(0))), ///
	vertical aseq swapnames ///
	coeflabels(L0 = "0" L1="1" L2="2" L3="3" L4="4" L5="5" L6="6") ///
	ylab(,nogrid) graphregion(color(white)) ///
	ytitle("Cumulative effect of 1-unit NAO increase" "on percentage of regional `outcome'") ///
	xtitle(Annual lags) graphregion(color(white)) legend(label(2 "`region'")) 
	graph export $figuresDir/fig2b-rev-`acro'.pdf, as(pdf) replace


//SOUTH ATLANTIC PLOT

	use $tempDir/select-sa, replace

	gen time1=year-1871
	gen time2 = time1^2
	gen time3=time1^3
	gen time4=time1^4
	gen time5=time1^5
	drop if year<1965

//FIG. 2A: SA CATCH PLOT
	gen l_catch=ln(main_mtons)
	
	local outcome "total catch"
	local acro "SA"
	local region "South Atlantic"
	local order "5"
	tsset year	

	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order', robust bw(auto)
	gen sample=1 if e(sample)
	ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm
	est store L0
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm
	est store L1
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm
	est store L2
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm
	est store L3
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm
	est store L4
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm
	est store L5
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm+L6.nao_slp_djfm
	est store L6
	
	coefplot ///
	(L0\L1\L2\L3\L4\L5\L6, recast(line) lcolor(dkgreen) lwidth(1) ciopts(recast(rarea) lcolor(green%0) fcolor(green%30) lwidth(1)) ///
	yline(0, lwidth(1) lc(gray) lp(dash)) ylabel(.15(.05)-.2, nogrid angle(0))), ///
	vertical aseq swapnames ///
	coeflabels(L0 = "0" L1="1" L2="2" L3="3" L4="4" L5="5" L6="6") ///
	ylab(,nogrid) graphregion(color(white)) ///
	ytitle("Cumulative effect of 1-unit NAO increase" "on percentage of regional `outcome'") ///
	xtitle(Annual lags) graphregion(color(white)) legend(label(2 "`region'")) 
	graph export $figuresDir/fig2a-catch-`acro'.pdf, as(pdf) replace

	drop sample _est_L0 _est_L1 _est_L2 _est_L3 _est_L4 _est_L5 _est_L6 l_catch

// FIG. 2B: SA REVENUE PLOT
	gen l_catch=ln(main_rev)
	local outcome "total revenue"

	tsset year
	
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order', robust bw(auto)
	gen sample=1 if e(sample)
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm
	est store L0
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)			
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm
	est store L1
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm
	est store L2
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm
	est store L3
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm
	est store L4
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm
	est store L5
	qui ivreg2  l_catch  L(0/6).nao_slp_djfm time1-time`order' if sample==1, robust bw(auto)
	lincomest L0.nao_slp_djfm+L1.nao_slp_djfm+L2.nao_slp_djfm+L3.nao_slp_djfm+L4.nao_slp_djfm+L5.nao_slp_djfm+L6.nao_slp_djfm
	est store L6

	coefplot ///
	(L0\L1\L2\L3\L4\L5\L6, recast(line) lcolor(dkgreen) lwidth(1) ciopts(recast(rarea) lcolor(green%0) fcolor(green%30) lwidth(1)) ///
	yline(0, lwidth(1) lc(gray) lp(dash)) ylabel(.15(.05)-.2, nogrid angle(0))), ///
	vertical aseq swapnames ///
	coeflabels(L0 = "0" L1="1" L2="2" L3="3" L4="4" L5="5" L6="6") ///
	ylab(,nogrid) graphregion(color(white)) ///
	ytitle("Cumulative effect of 1-unit NAO increase" "on percentage of regional `outcome'") ///
	xtitle(Annual lags) graphregion(color(white)) legend(label(2 "`region'")) 
	graph export $figuresDir/fig2b-rev-`acro'.pdf, as(pdf) replace


