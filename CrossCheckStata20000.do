clear
clear matrix
set mem 1g
set more off
insheet using CrossCheck2000.csv
xtset hhid tid
capture log close
log using CrossCheckStata.txt, replace
*timer on 1
*xi: xtreg y x1 x2 i.tid [aw=w], fe
*timer off 1
* Fixed effects, robust
timer on 2
xi: xtreg y x1 x2 i.tid [aw=w], fe vce(robust)
timer off 2
timer list
log close
