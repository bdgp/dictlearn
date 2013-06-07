function newStrVector = mergeTwoStrVectors(strVec1,strVec2,sep)
n = length(strVec1);
m = length(strVec2);
if n~=m
    sprintf('String vector lengths are not consistent!');
    newStrVector = [];
    return;
end
if ~exist('sep')
    sep = ' ';
end

newStrVector = cell(1,n);
for i = 1:n
    newStrVector{i} = [strVec1{i},sep,strVec2{i}];
end