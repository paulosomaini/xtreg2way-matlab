function s = regress1(y,X)
%REGRESS1 Regression
s.XX=X'*X;
s.beta=(s.XX)\(X'*y);
s.res=y-X*s.beta;
end

