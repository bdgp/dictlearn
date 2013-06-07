function Y = crossPat(width,height,crossW,crossH)
Y = zeros(height,width);
rowMidPt = floor(height/2);
colMidPt = floor(width/2);
rowInd = (rowMidPt-crossW):(rowMidPt + crossW);
colInd = (colMidPt+1):(colMidPt+crossH);
Y(rowInd,colInd) = 1;

rowInd = (rowMidPt+1):(rowMidPt + crossH);
colInd = (colMidPt-crossW):(colMidPt+crossW);
Y(rowInd,colInd) = 1;
