function [betaHat,aVarHat] = xtreg2wayPost(y,X,struc,se,noise)
%XTREG2WAY Estimates a 2-way fixed effect model absorbing the two set of
%dummies and reports standard errors. This function is meant to be used
%only after y and X have been projected on the set of dummies.
% Usage:
%  [betaHat,aVarHat] = xtreg2way(y,X,struc,se,noise)
%  Description of arguments
%  y (N-by-1) is the residual of the dependent variable after being
%  projected on the set of dummies.
%  X (N-by-K) is the matrix of residuals of covariates after being
%  projected on the set of dummies
%  struc (structure) contains results of the first step of the algorithm.
%  (obtained when calling xtreg2way).
%  se (possible values: 0,1,2,11) sets the standard error estimate. If
%  se==0 clasic standard errors assuming homoscedasticity and no within
%  group correlation or serial correlation. If se==1 standard errors
%  proposed by Arellano (1987) robust to heteroscedasticity and serial 
%  correlation. If se==2 it computes standard errors robust to heteroscedasticity, 
%  but assumes no correlation within group or serial correlation.
%  If se==11 Arellano (1987) standard errors with a degree of freedom
%  correction performed by Stata xtreg, fe. If se is omitted or set to [] 
%  then se is set to 1 and the Arellano (1987) estimator is computed.
%  noise (posible values 0,1) If noise==0, results are not displayed. If
%  noise==1 results are displayed. If noise is omitted or set to [],
%  results are displayed.
%  Description of output
%  betaHat (K-by-1) vector of estimated coefficients
%  aVarHat (K-by-K) estimate of the matrix of variances and covariance of 
%  the estimator.
%  Alternative Usage:
%  If the first step was already performed and variables yp and Xp are
%  the projected error of the original y and X on the two sets of
%  dummies, it is possible to call: 
% [betaHat,aVarHat] = xtreg2way(yp,Xp,struc);

if nargin<3, error('xtreg2way:nei','This function requires at least three inputs'); end
[obs,K]=size(X);
if obs~=numel(struc.hhid) || obs~=numel(struc.tid) || numel(y)~=obs
    error('xtreg2wayPost:dim',...
        'y, hhid and tid should be Nth order vectors, and X should be a N by K matrix')
end
if nargin<4 || isempty(se), se=1; end
if nargin<5 || isempty(noise), noise=1; end
reg=regress1(y,X);
betaHat=reg.beta';
switch se
    case 0
        sig2hat=(reg.res'*reg.res)/(sum(struc.w>0)-struc.N-struc.T+1-numel(reg.beta));
        aVarHat=sqrt(diag(sig2hat*inv(reg.XX)));
    case 1
        aVarHat=avar(X,reg.res,struc.hhid,reg.XX);
    case 2
        aVarHat=avar(X,reg.res,1:obs,reg.XX);
    case 11
        aVarHat=avar(X,reg.res,struc.hhid,reg.XX);
        stata_dof=((obs-1)/(obs-numel(reg.beta)-1))*(struc.N/(struc.N-1));
        aVarHat=aVarHat*(stata_dof)^2;
    otherwise
        disp('Computing standard errors robust to heteroskedasticity and within group correlation');
        aVarHat=avar(X,reg.res,struc.hhid,reg.XX);
end

if noise
    format('long');
    disp('Coefficient  S.E.     t-stat       p-val' )
    std=sqrt(diag(aVarHat));
    disp([betaHat'  std abs(betaHat'./std) (1-cdf('normal',abs(betaHat'./std),0,1))/2])
end
end

