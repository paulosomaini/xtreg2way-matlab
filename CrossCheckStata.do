clear
clear matrix
set mem 1g
set more off
insheet using CrossCheckBalanced1000.csv
xtset hhid tid
capture log close
log using CrossCheckStata.txt, replace
*XTREG
timer clear 
timer on 1
xi: xtreg y x1 x2 i.tid [aw=w], fe
timer off 1
*timer off 1
* Fixed effects, robust
timer on 2
xi: xtreg y x1 x2 i.tid [aw=w], fe vce(robust)
timer off 2
*REGHDFE
timer on 3
reghdfe y x1 x2, absorb(hhid tid)
timer off 3
* Fixed effects, robust
timer on 4
reghdfe y x1 x2, absorb(hhid tid) vce(cluster hhid)
timer off 4
timer list
log close