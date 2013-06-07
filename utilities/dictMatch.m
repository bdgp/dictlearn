function matchInd = dictMatch(Dtemplate,D,method)
% Function that matches a new dict with a template
% Dtemplate: a template to be aligned with
% D: the incoming dictionary
% method: distance (or alike) metric
disMat = dissimilarityDict(Dtemplate,D,method);
[matchInd,cost] = munkres(disMat);