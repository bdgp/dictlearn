function Dstd = dictStd(D,method)
% D: dictionary
% method: 1: by L-1 norm
% 2: by L-2 norm
% 0: by max
[p,K] = size(D);
Dstd = D;
for i = 1:K
    if method == 0
        Dstd(:,i) = D(:,i)/max(D(:,i));
    elseif method == 1
        l1Norm = sum(abs(D(:,i)));
        Dstd(:,i) = D(:,i)/l1Norm;
    elseif method == 2
        l2Norm = norm(D(:,i));
        Dstd(:,i) = D(:,i)/l2Norm;
    end
end

