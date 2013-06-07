function disMatrix = dissimilarityDict(D1,D2,method)
% This function computes pairwise dissimilarity measure of
% dictionary elements in D1 and D2
% method: 'corr': the 1 - correlation measure
%         'cos': cosine of the pairwise angles
%         'euclidean': square euclidean distance
addpath('~/fruitfly/annotation/TF/annotSOM/somtoolbox/pwmetric/');
if size(D1,1)~=size(D2,1) | size(D1,2)~=size(D2,2)
    sprintf('The sizes of the two dictionaries are not consistent!');  
    disMatrix = [];
    return
end

if strmatch(method,'corr','exact')
    cor = corr(D1,D2);
    disMatrix = 1 - cor;    
elseif strmatch(method,'cos','exact')
    disMatrix = slmetric_pw(D1,D2,'corrdist');
elseif strmatch(method,'euclidean')
    disMatrix = slmetric_pw(D1,D2,'sqdist');
elseif strmatch(method,'regression')
    param.lambda = 0;
    param.mode = 2;
    param.pos = true;
    A = mexLasso(D1,D2,param);
    B = mexLasso(D2,D1,param);
    disMatrix = amariError(A) + amariError(B);
end

