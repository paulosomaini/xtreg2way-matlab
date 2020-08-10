function [betaHat,aVarHat,y,X,struc,cluster] = xtreg2way(y,X,iid,tid,w,struc,se,cluster,noise)
%XTREG2WAY Estimates a 2-way fixed effect model absorbing the two set of
%dummies and reports standard errors
% Usage:
%  [betaHat,aVarHat,y,X,struc] = xtreg2way(y,X,iid,tid,w,struc,se,noise)
%  Description of arguments
%  y (N-by-1) is the dependent variable
%  X (N-by-K) is the matrix of covariates
%  iid (N-by-1) is the group id
%  tid (N-by-1) is the time id
%  w (N-by-1) weights. If w is omitted or set to [] then w=1 for all obs
%  struc (structure) contains the results of the first step of the
%  algorithm. If struc is omitted or set to [] then the first step is
%  performed.
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
%  y (N-by-1) the residual of the projection of y on the two sets of
%  dummies.
%  X (N-by-K) the residual of the projection of each column of X on the two
%  sets of dummies.
%  struc (structure) results of the first step of the algorithm.
%  Alternative Usage:
%  If the first step was already performed and variables yp and Xp are
%  the projected error of the original y and X on the two sets of
%  dummies, it is possible to call: 
% [betaHat,aVarHat] = xtreg2way(yp,Xp,struc);

if nargin<3, error('xtreg2way:nei','This function requires at least three inputs'); end
[obs,K]=size(X); projectVars=1;
if nargin<4, projectVars=0; struc=iid; iid=struc.hhid; tid=struc.tid; w=struc.w; end
if obs~=numel(iid) || obs~=numel(tid) || numel(y)~=obs
    error('xtreg2way:dim',...
        'y, iid and tid should be Nth order vectors, and X should be a N by K matrix')
end
if nargin<5 || isempty(w), w=ones(obs,1); end
[flag_redundant,nr] = nonredundant(iid,tid,w);
if flag_redundant
    esample = ismember(iid,nr.iid) & ismember(tid,nr.tid);
    y = y(esample); X = X(esample,:); iid = iid(esample); 
    tid = tid(esample);w = w(esample);
end
if nargin<6 || isempty(struc), struc=projdummies(iid,tid,w); end
if nargin<7 || isempty(se), se=1; end
if nargin<8 || isempty(cluster), cluster=struc.hhid; end
if nargin<9 || isempty(noise), noise=1; end

if flag_redundant, struc.esample = esample; end
if projectVars
    for kk=1:K,X(:,kk)=projvar(X(:,kk),struc);end
    y=projvar(y,struc);
end
reg=regress1(y,X);
betaHat=reg.beta';
dof =struc.obs /(struc.obs-struc.N-struc.T-numel(reg.beta));
switch se
    case 0
        sig2hat=(reg.res'*reg.res)/(sum(struc.w>0)-struc.N-struc.T+1-numel(reg.beta));
        aVarHat=sig2hat*inv(reg.XX);
    case 1
        aVarHat=avar(X,reg.res,cluster(esample,:),reg.XX)*dof;
    case 2
        aVarHat=avar(X,reg.res,1:obs,reg.XX)*dof;
    case 11
        aVarHat=avar(X,reg.res,cluster(esample,:),reg.XX);
        stata_dof=((obs-1)/(obs-numel(reg.beta)-1))*(struc.N/(struc.N-1));
        aVarHat=aVarHat*(stata_dof)^2;
    otherwise
        disp('Computing standard errors robust to heteroskedasticity and within group correlation');
        aVarHat=avar(X,reg.res,struc.hhid,reg.XX)*dof;
end

if noise
    format('short');
    disp('Coefficient  S.E.     t-stat       p-val' )
    std=sqrt(diag(aVarHat));
    disp([betaHat'  std abs(betaHat'./std) (1-cdf('normal',abs(betaHat'./std),0,1))/2])
end
end

