## An Algorithm to Estimate the Two-Way Fixed Effect Model

This folder contains the Matlab functions necessary to run the algorithm described in the notes. It also contains a matlab script that generates a dataset and runs the code, and a Stata do file that runs the same procedure on the same data using xtreg and adding time dummies.

Wrapper:

`reg2way.m` estimates a 2-way fixed effect model absorbing the two set of dummies and reports standard errors. 

Syntax:

```
[betaHat,aVarHat,yp,Xp,struc,cluster]=xtreg2way(y,X,iid,tid,w,struc,se,cluster,noise)
```
## Description of arguments

- `y` (N-by-1) is the dependent variable
- `X` (N-by-K) is the matrix of covariates
- `iid` (N-by-1) is the group id
- `tid` (N-by-1) is the time id
- `w` (N-by-1) is a vector of weights. If w is omitted or set to [] then w=1 for all observations.
- `struc` (structure) contains the results of the first step of the  algorithm. struc can be omitted or set to [] in which case, it will be automatically calculated. Notice that the first step is computationally expensive. If many regressions are to be run using the same set of covariates it is possible to run the first step only once. 
- `se` (possible values: 0,1,2) 
	- `se`==0 : standard errors assuming homoscedasticity and no within  group correlation or serial correlation. This option can not be used with cluster option.
	- `se`==1 : standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation. This option works with cluster option. If cluster option is omitted, the the standard errors will be calculated clustering by the first fixed effect selected. 
	- `se`==2 : standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation. This option can not be used with cluster option.
- `cluster` variable with cluster id. If `se` is omitted or set to [] and `cluster` is selected, then `se` is set to `1`. If `se` and `cluster` are omitted or set to [], then `se` is set to `"1"` and `cluster` is set to `iid`, and the Arellano (1987) estimator is computed.  
-`noise` (possible values 0,1):
	- `noise`==0 results are not displayed. 
	- `noise`==1 results are displayed (default). 

## Description of output:

- `betaHat` (K-by-1) vector of estimated coefficients
- `aVarHat` (K-by-K) estimate of the matrix of variances and covariance of  the estimator.
- `yp` (N-by-1) the residual of the projection of y on the two sets of  dummies.
- `Xp` (N-by-K) the residual of the projection of each column of X on the two  sets of dummies.
- `struc` (structure) results of the first step of the algorithm.
- `cluster` make an esample of the variable of the cluster selected. 

## Alternative Usage:

If the first step was already performed and variables `yp` and `Xp` are  the projected error of the original `y` and `X` on the two sets of  dummies, then it is possible to avoid repeating the first two steps of the algorithm by calling:

```[betaHat,aVarHat]=xtreg2way(yp,Xp,struc)  ```

or

```[betaHat,aVarHat]=xtreg2wayPost(yp,Xp,struc,se,cluster,noise)```

Matlab Functions:

- projdummies.m computes the first step of the algorithm. Let S=[D,H], where D is the matrix of individual effect dummies and H is the matrix of time dummies. It does not compute the inverse of (S’S), but it returns a structure with the minimal information required to construct it. Syntax: `struc=projdummies(hhid,tid,w)`. Observations with weights equal to zero are dropped. The three vectors have to have the same length. 
- projvar.m obtains the residualized variables. It uses the structure generated by projdummies. Syntax: `yp=projvar(y,struc)`, where `yp` is the residualized version of `y`. Also: `Xp=projvar(X,struc)`, where `X` and `Xp` are matrices.
- regress1.m performs an OLS regression (third step) and stores a structure with coefficients, residuals and the matrix `Xp'Xp`. Syntax: `reg=regress1(yp,Xp)`, where `reg` is a structure, `yp` is the result of `projvar(y,struc)` and `Xp` is the result of stacking horizontally all vectors `projvar(x,struc)` for each `x` that is a column of `X`. The structure `reg` has three fields: `beta`, `res`, and `XX` which respectively stand for coefficients, residuals, and the matrix `(Xp'Xp)`.
- avar.m computes the asymptotic variance of the estimator (fourth step). Syntax: `matCov=avar(Xp,e,group,J)`. 
	- `Xp` is an L-by-K matrix of residualized covariates, where L is the number of observations, and K is the number of covariates. 
	- `e` is the L-by-1 vector of residuals. 
	- `group` is the L-by-1 cluster identifier (it should coincide with hhid in step one if standard errors are clustered at the individual level. If `group` has more or less than L elements it will be assumed that the user required the heteroscedasticity robust estimate. Otherwise, the clustered heteroscedasticity robust estimate is provided. 
	- If the matrix `Xp'Xp` was calculated before, it can be entered as the fourth argument. `J` has to be a K-by-K matrix.

## References

Arellano, M. (1987), Computing Robust Standard Errors for Within-Groups Estimators, Oxford Bulletin of Economics and
Statistics, 49, issue 4, p. 431–434. 

Cameron, A. C., & Miller, D. L. (2015). A Practitioner’s Guide to Cluster-Robust Inference. Journal of Human Resources, 50(2), 317-372.

Somaini, P. and F.A. Wolak, (2015), An Algorithm to Estimate the Two-Way Fixed Effects Model, Journal of Econometric Methods, 5, issue 1, p. 143-152.


