function strVector = num2strBatch(numVector)
n = length(numVector);
strVector = cell(1,n);
for i = 1:n
    temp = round(100*numVector(i))/100;
    strVector{i} = num2str(temp);
end
    