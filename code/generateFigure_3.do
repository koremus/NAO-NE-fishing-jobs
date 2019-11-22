/*******************************************************************************
generateFigure_3
Script Description: This STATA script inputs and merges data for regression. 
Code developed for STATA v15.

Author:  Kimberly L. Oremus
Date: April 29, 2019

Data from: 
Bureau of Labor Statistics
NAICS 11411 Fishing (and all its subindustries 114111, 114112, 114119)
42446 FishSeafoodWholesalers
44522 FishSeafoodMarkets
52 FinanceInsurance
61 Education
71 Arts and Recreation
5615 Travel
56152 TourOperators
21 Extraction
Many other industries were also looked at, but do not appear in the article due to space constraints.		
Documentation on NAICS: https://www.census.gov/eos/www/naics/2012NAICS/2012_Definition_File.pdf
Due to the size of the bls files, the data has already been cleaned and processed and only includes
the following coastal states: 
New England (NE): ME, NH, MA, RI, CT	
Mid Atlantic (MA): VA, NY, MD, PA, NJ, DE 
South Atlantic (SA): SC, GA, FL, NC
NAO data (described below) was merged in a similar fashion to the script in generateFigure_1_2.do
data located in input folder "bls.dta"
NAICS and SIC data downloaded March 12, 2019


NCAR ClimateDataGuide
DJFM annual NAO index
https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based
data located in input folder "NAO_slp_djfm_1864_2017.csv"
accessed on January 10, 2019
*******************************************************************************/

clear all
set trace off
set more off
capture log close

cap ssc install coefplot
cap ssc install lincomest
cap ssc install estout
cap ssc install ranktest 
cap ssc install reghdfe
cap ssc install ftools
cap ssc install ivreg2 


//set local directory 
cd ..

global inputDir "input"
global resultsDir "tables"
global figuresDir "figures"
global tempDir "temp"

use $inputDir/bls, replace
drop if MA==1

//totalW denotes total wages
//emp donotes employment 
//estabs denotes number of establishments
foreach f in totalW emp estabs {

gen var=`f'Fishing

// turn data into hyperbolic sine
gen fish=log(var + sqrt(var^2 + 1))

// generate balanced sample
bysort fipscode: gen temp = _N if fish != . 
bysort fipscode: egen obsN = max(temp)
drop temp 

summ obsN 
gen balancedSample = (obsN == r(max)) 
gen treatment=NE*nao_slp_djfm

xtset fipscode year

******************* Fig. 3: Fishing Employment ************************/

qui reghdfe fish L(0/6).treatment if balancedSample==1, a(fipscode year) cluster(fipscode)
gen sample=1 if e(sample)
qui eststo M1: reghdfe fish L(0/6).treatment if sample==1, a(fipscode year) cluster(fipscode)
drop sample

gen econ=log(`f'Travel + sqrt(`f'Travel^2 + 1))
bysort fipscode: gen temp = _N if econ != . 
bysort fipscode: egen obsN2 = max(temp)
drop temp 
summ obsN2 
gen balancedSample2 = (obsN2 == r(max)) 
qui reghdfe fish L(0/6).treatment if balancedSample==1 & balancedSample2==1, a(fipscode year) cluster(fipscode)
gen sample=1 if e(sample)
qui eststo M2: reghdfe econ L(0/6).treatment if sample==1, a(fipscode year) cluster(fipscode)
drop obsN2 econ balancedSample2 sample

gen econ=log(`f'FinanceInsurance + sqrt(`f'FinanceInsurance^2 + 1))
bysort fipscode: gen temp = _N if econ != . 
bysort fipscode: egen obsN2 = max(temp)
drop temp 
summ obsN2 
gen balancedSample2 = (obsN2 == r(max)) 
qui reghdfe fish L(0/6).treatment if balancedSample==1 & balancedSample2==1, a(fipscode year) cluster(fipscode)
gen sample=1 if e(sample)
qui eststo M3: reghdfe econ L(0/6).treatment if sample==1, a(fipscode year) cluster(fipscode)
drop obsN2 econ balancedSample2 sample

gen econ=log(`f'Education + sqrt(`f'Education^2 + 1))
bysort fipscode: gen temp = _N if econ != . 
bysort fipscode: egen obsN2 = max(temp)
drop temp 
summ obsN2 
gen balancedSample2 = (obsN2 == r(max)) 
qui reghdfe fish L(0/6).treatment if balancedSample==1 & balancedSample2==1, a(fipscode year) cluster(fipscode)
gen sample=1 if e(sample)
qui eststo M4: reghdfe econ L(0/6).treatment if sample==1, a(fipscode year) cluster(fipscode)
drop obsN2 econ balancedSample2 sample

gen econ=log(`f'ArtsRecreation + sqrt(`f'ArtsRecreation^2 + 1))
bysort fipscode: gen temp = _N if econ != . 
bysort fipscode: egen obsN2 = max(temp)
drop temp 
summ obsN2 
gen balancedSample2 = (obsN2 == r(max)) 
qui reghdfe fish L(0/6).treatment if balancedSample==1 & balancedSample2==1, a(fipscode year) cluster(fipscode)
gen sample=1 if e(sample)
qui eststo M5: reghdfe econ L(0/6).treatment if sample==1, a(fipscode year) cluster(fipscode)
drop obsN2 econ balancedSample2 sample

gen econ=log(`f'Extraction + sqrt(`f'Extraction^2 + 1))
bysort fipscode: gen temp = _N if econ != . 
bysort fipscode: egen obsN2 = max(temp)
drop temp 
summ obsN2 
gen balancedSample2 = (obsN2 == r(max)) 
qui reghdfe fish L(0/6).treatment if balancedSample==1 & balancedSample2==1, a(fipscode year) cluster(fipscode)
gen sample=1 if e(sample)
qui eststo M6: reghdfe econ L(0/6).treatment if sample==1, a(fipscode year) cluster(fipscode)
 
coefplot ///
(M2, keep(*treatment) mcolor(gs10%0) ciopts(lcolor(gs10%0) lwidth(.5))) ///
(M2, keep(*treatment) mcolor(gs12) ciopts(lcolor(gs12) lwidth(.5))) ///
(M3, keep(*treatment) mcolor(gs10) ciopts(lcolor(gs10) lwidth(.5))) ///
(M4, keep(*treatment) mcolor(gs8) ciopts(lcolor(gs8) lwidth(.5))) ///
(M5, keep(*treatment) mcolor(gs6) ciopts(lcolor(gs6) lwidth(.5))) ///
(M6, keep(*treatment) mcolor(purple) ciopts(lcolor(purple) lwidth(.5))) ///
(M1, keep(*treatment) mcolor(eltblue) ciopts(lcolor(eltblue) lwidth(.5)) ///
yline(0, lwidth(1) lc(black)) ///
graphregion(color(white)) ///
ylabel(, nogrid angle(0))) ///
(M2, keep(*treatment) mcolor(gs10%0) ciopts(lcolor(gs10%0) lwidth(.5))) ///
(M2, keep(*treatment) mcolor(gs10%0) ciopts(lcolor(gs10%0) lwidth(.5))), ///
xtitle("Annual lags") ytitle("% Effect of 1-unit NAO increase on `f'") ///
vertical coeflabels(treatment="0" L.treatment="1" L2.treatment="2" L3.treatment="3" ///
L4.treatment="4" L5.treatment="5" L6.treatment="6") ylab(,nogrid) graphregion(color(white)) ///
legend(off) 
graph export $figuresDir/fig3_`f'fishing_fe.pdf, as(pdf) replace
set trace off
set more off
capture log close

drop econ obsN obsN2 balancedSample balancedSample2 treatment sample var fish
}

