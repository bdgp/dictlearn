function sqmsk = squareMask(width,height,blockSize)
sqmsk = zeros(height,width);
ncol = width/blockSize;
nrow = height/blockSize;
rowon = 1;
for i = 1:nrow
    if rowon == 1
        colon = 1;
    else
        colon = 0;
    end
    for j = 1:ncol
        if colon
            rowInd = (1+blockSize*(i-1)):(blockSize*i);
            colInd = (1+blockSize*(j-1)):(blockSize*j);
            sqmsk(rowInd,colInd) = 1;
        end
        colon = 1 - colon;
    end
    rowon = 1 - rowon;
end