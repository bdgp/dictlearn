function D = eroImg(d,level,width,height)
D = zeros(size(d,1),level);
temp = reshape(d,height,width);

se = strel('line',5,0);
for i = 1:level
    temp = imerode(temp,se);
    D(:,i) = temp(:);
end
