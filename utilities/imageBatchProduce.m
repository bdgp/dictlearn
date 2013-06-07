function I = imageBatchProduce(D,width,height,nrow,ncol)
% D: a list of images... each column corresponding to an image
% width/height: the width/height of an image
% nrow/ncol: how many rows/cols you want to have in the layout?
% colorScale: can be '' or 'gray'
    [p,n] = size(D);
   nR = nrow*height;
   nC = ncol*width;
   I = zeros(nR,nC);
   rowStart = 1;
   colStart = 1;
   k = 1;
   for i = 1:nrow
       for j = 1:ncol
           indR = (rowStart + (i-1)*height):( i*height);
           indC = (colStart + (j-1)*width):( j*width);
           I(indR,indC) = reshape(D(:,k),height,width);
           k = k+1;
           if k>n
               break;
           end
       end
       if k>n
           break;
       end
   end
end
