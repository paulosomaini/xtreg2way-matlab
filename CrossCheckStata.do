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
timer list 
timer clear 
outreg2 using XTRegTestSTATA1000.doc, replace ctitle(Xtreg) keep(x1 x2) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1))
* Fixed effects, robust
timer on 1
xi: xtreg y x1 x2 i.tid [aw=w], fe vce(robust)
timer off 1
timer list 
timer clear 
outreg2 using XTRegTestSTATA1000.doc, append ctitle(Xtreg-R) keep(x1 x2) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1)) 
timer clear 
*Change specification
timer on 1
xi: xtreg y x1 i.tid [aw=w], fe
timer off 1
timer list 
timer clear 
outreg2 using XTRegTestSTATA1000.doc, append ctitle(Xtreg-one C) keep(x1) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1))
timer on 1
xi: xtreg y x1 i.tid [aw=w], fe vce(robust)
timer off 1
timer list 
timer clear 
outreg2 using XTRegTestSTATA1000.doc, append ctitle(Xtreg-one C-R) keep(x1) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1))

*REGHDFE
timer on 1
reghdfe y x1 x2, absorb(hhid tid)
timer off 1
timer list 
timer clear 
outreg2 using ReghdfeTestSTATA1000.doc, replace ctitle(Reghdfe) keep(x1 x2) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1))
* Fixed effects, robust
timer on 1
reghdfe y x1 x2, absorb(hhid tid) vce(cluster hhid)
timer off 1
timer list 
timer clear 
outreg2 using ReghdfeTestSTATA1000.doc, append ctitle(Reghdfe-R) keep(x1 x2) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1))
*change specification
timer on 1
reghdfe y x1, absorb(hhid tid)
timer off 1
timer list 
timer clear 
outreg2 using ReghdfeTestSTATA1000.doc, append ctitle(Reghdfe-one C) keep(x1) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1))
* Fixed effects, robust
timer on 1
reghdfe y x1, absorb(hhid tid) vce(cluster hhid)
timer off 1
timer list 
timer clear 
outreg2 using ReghdfeTestSTATA1000.doc, append ctitle(Reghdfe-one C-R) keep(x1) addtext(Group FE, YES, Time FE, YES) addstat(Running Time, r(t1))

log close