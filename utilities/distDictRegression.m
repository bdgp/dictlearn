function [dist1, dist2, error1, error2] = distDictRegression(D1,D2)
% This function gives a quantitative measure of two dictionaries
% the items in both dictionary may not be in the same column and
% thus a direct calculation is not feasible.
% This function uses instead a regression approach, by regressing
% D1 on D2 first to get the regression coefficient matrix. If two
% dictionaries are similar the regression coefficient matrix should
% be close to an identity matrix after a proper permutation of the 
% columns. The closeness to an identity matrix is measured by Amari 
% metric, which is commonly used in ICA.
% error1, error2: relative representation error.
param.lambda = 0;
param.mode = 2;
param.pos = true;
A = mexLasso(D1,D2,param);
B = mexLasso(D2,D1,param);

dist1 = amariError(A);
E1 = D1 - D2*A;
error1 = mean(sum(E1.^2)./sum(D1.^2));

dist2 = amariError(B);
E2 = D2 - D1*B;
error2 = mean(sum(E2.^2)./sum(D2.^2));

