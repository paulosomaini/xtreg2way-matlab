function matCov=avar(X,e,group,J)
%FUNCTION avar(X,e,group,J)
%Computes the asymptotic variance of the OLS estimator. 
%X is the matrix of covariates: L-by-K. L num of obs, and K num of
%covariates.
%e is the vector of residuals: L-by-1;
%group is the cluster level: L-by-1. If g has more or less than L elements
%it will be assumed that the user required the heteroscedasticity robust
%estimate. Otherwise, the clustered heteroscedasticity robust estimate is
%provided.
% If the matrix X'*X was calculated before, it can be entered as the
% argument J. J is a K-by-K matrix.
[L,K]=size(X); e=e(:); L2=numel(e);
if ~(L==L2),error('myApp:dimen','X and e should have the same number of rows'); end
if nargin<4;
    J=X'*X;
end
[~,~,group]=unique(group); G=max(group);
eX=sparse(1:L,1:L,e,L,L,L)*X;
if nargin<3 || ~(numel(group)==L) || G==L;
    disp('Heteroscedasticity Robust SE');
    V=eX'*eX;
else
    eX=sparse(group,1:L,1,G,L,L)*eX;
    V=eX'*eX;
% Previous slower version:
%     V2=zeros(K,K);    
%     for g=1:G
%         Xi=X(group==g,:);
%         ei=e(group==g,:);
%         Xe=Xi'*ei;
%         V2=V2+Xe*Xe';
%     end
%     disp(max(abs(V(:)-V2(:))));
%     disp(mean(V(:)));

end
matCov=J\V/J;
