**  DO-FILE METADATA
//  algorithm name				scdgrowth_001.do
//  project:							SCD growth and effect on later life outcomes
//  analysts:							Ian HAMBLETON, Christina HOWITT
//	date last modified		11-Aug-2018
//  algorithm task				Load cohort data. Generate growth z-scores

** General algorithm set-up
version 15
clear all
macro drop _all
set more 1
set linesize 80

** Set working directories: this is for DATASET and LOGFILE import and export
** DATASETS to encrypted SharePoint folder
local datapath "X:\The University of the West Indies\DataGroup - repo_data\data_p123"
** LOGFILES to unencrypted OneDrive folder
local logpath X:\OneDrive - The University of the West Indies\repo_datagroup\repo_p123

** Close any open log fileand open a new log file
capture log close
cap log using "`logpath'\scdgrowth_001", replace

** In this project, we:
** 	(1) Explore whether SS children are wasted compared to WHO growth standards.
**	(2) We then use these standardised growth metrics to explore their effect
**			on a series of later-life health outcomes
tempfile one two three_ss three_aa four_ss four_aa five

** (1) 	Loading Static Information from Jamaican Cohort Study
**			We use
**			genotype 						--> geno
**			cohort indicator 		--> cohort
**			ID Number						--> idno
**			Date of last visit	--> dolv
** use "X:\The University of the West Indies\DataGroup - repo_data\data_p123\all_static", clear
use "`datapath'\all_static", clear
keep idno cohort sex geno status dob dolv

** Genotype restriction, SS==1 AA==2
recode geno 29=1
keep if cohort==1 & geno<=2

** Sex (for zanthro command) must be 1 and 2
recode sex 0=1 1=2
label define sex_ 1 "male" 2 "female",modify
label values sex sex_

** Drop 4 SS patients default without ever attending clinic
drop if idno==1 | idno==6 | idno==13 | idno==3163

** DECISION. Drop participants without a minimum of 2 years of follow-up
** This is somewhat arbitrary when using the WHO external standard
** Which allows measurement from birth...
** How many patients are lost to follow up before 2 years of age?
gen ay = (dolv-dob)/365.25
mark out if ay<2
egen tag = tag(idno)
count if tag==1
tab out status if tag==1 & geno==1
tab out status if tag==1 & geno==2
label var out "Less than 2 yrs of f/u"
label var ay "Age at last visit in years"
** Save temporary file
sort idno
save `one'

** Merge with Jamaican cohort height and weight information
use "`datapath'\all_hw", clear
sort idno dov
merge idno using `one'
keep if _merge==3
drop _merge tag
sort idno dov
save `two'

** Generate age in years and months
gen agey = (dov-dob)/365.25
gen agem = agey*12
label var agey "Age in years at each measurement"
label var agem "Age in months at each measurement"

** Reporting on childhood and adolescence (2-18 years)
** We also restrict out measurements before 2 years
** This is the minimum age for international BMI standards
** 25 months = 2 years 1 months and up
** 228 months = 19 years
keep if agem>=25 & agem<228

** For graphics we use AA and SS
** For outcome modelling we will restrict to SS only
keep if geno<=2
label define geno 1 "ss" 2 "aa",modify
label values geno geno

** Create 6-month age groups for comparison with
** external reference population
** 24 months    = 2 years
** 216.5 months = 18 years
gen agec = recode(agem,24,24.5,30.5,36.5,42.5,48.5,54.5,60.5,66.5,  /*
*/                     72.5,78.5,84.5,90.5,96.5,                    /*
*/                     102.5,108.5,114.5,120.5,126.5,132.5,         /*
*/                     138.5,144.5,150.5,156.5,162.5,168.5,         /*
*/                     174.5,180.5,186.5,192.5,198.5,204.5,         /*
*/                     210.5,216.5)

** Drop unecessary variables
keep idno sex agey agem agec ht wt geno
order idno sex geno agey agec agem ht wt

** Choosing nearest ht and wt measurement to year mid-point
gen d1 = agem-agec
gen d2 = (d1^2)^(1/2)
egen d3 = min(d2), by(idno agec)
keep if d3==d2

** Ensure just 1 measurement per person per year
** n should contain only 1's
sort idno agec
by idno agec: gen n = _n

** Calculate BMI-for-age in JSSCD sample
gen bmi = wt / ((ht/100)^2)

** ----------------------------------------------------
** 11-AUG-2018
** Drop POTENTIAL ERRORS - Height and Weight measurement.
** This is a temporary clean for the proposal deadline
** Real cleaning needed
** ----------------------------------------------------
drop if bmi>30 | bmi<11

** Use EGEN to create external reference centiles
** SS only first
preserve
 	keep if geno==1
	** BMI-for-age
	egen zbmi_who= zanthro(bmi, ba, WHO), xvar(agem) gender(sex) gencode(m=1, f=2) ageunit(month)
	egen zbmi_us= zanthro(bmi, ba, US), xvar(agem) gender(sex) gencode(m=1, f=2) ageunit(month)
	egen zbmi_uk= zanthro(bmi, ba, UK), xvar(agem) gender(sex) gencode(m=1, f=2) ageunit(month)
	** Save this file
	save `three_ss'

	** Average BMI across participants
	collapse (mean) bmi, by(geno sex agec)
	save `four_ss'
restore

** AA only
preserve
	keep if geno==2
	** BMI-for-age
	egen zbmi_who= zanthro(bmi, ba, WHO), xvar(agem) gender(sex) gencode(m=1, f=2) ageunit(month)
	egen zbmi_us= zanthro(bmi, ba, US), xvar(agem) gender(sex) gencode(m=1, f=2) ageunit(month)
	egen zbmi_uk= zanthro(bmi, ba, UK), xvar(agem) gender(sex) gencode(m=1, f=2) ageunit(month)
	** Save this file
	save `three_aa'

	** Average BMI across participants
	collapse (mean) bmi, by(geno sex agec)
	save `four_aa'
restore

** Load WHO standard for girls
use "`datapath'\who_zscore_girls", clear
gen sex = 2
append using "`datapath'\who_zscore_boys"
replace sex = 1 if sex==.
replace month = month+0.5
rename month agec
save `five'

** Load the SCD datasets and merge  the WHO standard without standardizing variable names
** We will only use this to create background curves
use `four_ss', clear
append using `four_aa'
** WHO standard data
merge m:1 sex agec using `five'


** GRAPHIC
** MEN -- BMI-for-age
#delimit ;
graph twoway  (rarea sd2neg sd1 agec if sex==1, sort fcolor("26 152 80%30") lw(none) )
							(rarea sd3neg sd2neg agec if sex==1, sort fcolor("215 48 39%30") lw(none) )
							(rarea sd1 sd2 agec if sex==1, sort fcolor("215 48 39%30") lw(none) )
							(line sd0 agec if sex==1, sort clw(0.25) clp("l") clc(gs12))

							(line bmi agec if agec>=60 & geno==1 & sex==1, sort clw(0.25) clp("l") clc(gs0))
							(line bmi agec if agec>=60 & geno==2 & sex==1, sort clw(0.25) clp("-") clc(gs0))
							,
							plotregion(c(gs16) ic(gs16) ilw(thin) lw(thin))
							graphregion(color(gs16) ic(gs16) ilw(thin) lw(thin))
							ysize(20) xsize(15)

       				xlab(36 "3" 72 "6" 108 "9" 144 "12" 180 "15" 216 "18", labs(3) nogrid glc(gs14))
       				xtitle("Age (years)", margin(t=3) size(3))
       				xscale(lw(vthin) range(24 300))

       				ylab(, labs(3) nogrid glc(gs14) angle(0) format(%9.0f))
							///ytick(10(20)80)
       				ytitle("BMI", margin(r=3) size(3))
       				yscale(lw(vthin))
							ytick(10(1)30)

							text(16.5 235 "Thinness", place(e) size(small) c(gs6))
							text(14.5 235 "Severe thinness", place(e) size(small) c(gs6))
							text(21.0 235 "Normal", place(e) size(small) c(gs6))
							text(27.0 235 "Overweight", place(e) size(small) c(gs6))
							///text(67.5 219 "50th", place(e) size(small) c(gs9))
							///text(75.5 219 "75th", place(e) size(small) c(gs9))

							legend(bexpand size(3) position(12) bm(t=1 b=0 l=0 r=0) colf cols(1)
							region(fcolor(gs16) lw(vthin) margin(l=2 r=2 t=2 b=2) )
							order(5 6)
							lab(5 "SS participants")
							lab(6 "AA participants")
							)
							name(men)
							;
#delimit cr


** GRAPHIC
** WOMEN -- BMI-for-age
#delimit ;
graph twoway  (rarea sd2neg sd1 agec if sex==2, sort fcolor("26 152 80%30") lw(none) )
							(rarea sd3neg sd2neg agec if sex==2, sort fcolor("215 48 39%30") lw(none) )
							(rarea sd1 sd2 agec if sex==2, sort fcolor("215 48 39%30") lw(none) )
							(line sd0 agec if sex==2, sort clw(0.25) clp("l") clc(gs12))

							(line bmi agec if agec>=60 & geno==1 & sex==2, sort clw(0.25) clp("l") clc(gs0))
							(line bmi agec if agec>=60 & geno==2 & sex==2, sort clw(0.25) clp("-") clc(gs0))
							,
							plotregion(c(gs16) ic(gs16) ilw(thin) lw(thin))
							graphregion(color(gs16) ic(gs16) ilw(thin) lw(thin))
							ysize(20) xsize(15)

       				xlab(36 "3" 72 "6" 108 "9" 144 "12" 180 "15" 216 "18", labs(3) nogrid glc(gs14))
       				xtitle("Age (years)", margin(t=3) size(3))
       				xscale(lw(vthin) range(24 300))

       				ylab(, labs(3) nogrid glc(gs14) angle(0) format(%9.0f))
							///ytick(10(20)80)
       				ytitle("BMI", margin(r=3) size(3))
       				yscale(lw(vthin))
							ytick(10(1)30)

							text(15.5 235 "Thinness", place(e) size(small) c(gs6))
							text(13.5 235 "Severe thinness", place(e) size(small) c(gs6))
							text(21.0 235 "Normal", place(e) size(small) c(gs6))
							text(27.0 235 "Overweight", place(e) size(small) c(gs6))
							///text(67.5 219 "50th", place(e) size(small) c(gs9))
							///text(75.5 219 "75th", place(e) size(small) c(gs9))

							legend(bexpand size(3) position(12) bm(t=1 b=0 l=0 r=0) colf cols(1)
							region(fcolor(gs16) lw(vthin) margin(l=2 r=2 t=2 b=2) )
							order(5 6)
							lab(5 "SS participants")
							lab(6 "AA participants")
							)
							name(women)
							;
#delimit cr

** Save final dataset
use `three_ss', replace
drop d1 d2 d3 n
label var agec "6-month age categories"
label var bmi "Body mass index"
label var zbmi_who "WHO bmi-for-age standard"
label var zbmi_us "US bmi-for-age standard"
label var zbmi_uk "UK bmi-for-age standard"
label data "Jamaican Cohort childhood BMI"
save "`datapath'\scdgrowth_001", replace
