function [distM,rowTemp,colTemp] = amariError(A)
% amari error
% normalized to [0,1]

[n,m] = size(A);
if n == 1 & m == 1
    distM = double(A==0);
    rowTemp = 0;
    cowTemp = 0;
    return;
end

maxCol = max(abs(A),[],1);
sumCol = sum(abs(A),1);
colTemp0 = abs(sumCol./maxCol - 1);
naInd = find(isnan(colTemp0)==1);
colTemp0(naInd) = n - 1;
colTemp = mean(colTemp0)/(n-1);


maxRow = max(abs(A),[],2);
sumRow = sum(abs(A),2);
rowTemp0 = abs(sumRow./maxRow - 1);
naInd = find(isnan(rowTemp0)==1);
rowTemp0(naInd) = m - 1;
rowTemp = mean(rowTemp0)/(m-1);

distM = (rowTemp + colTemp)/(2);