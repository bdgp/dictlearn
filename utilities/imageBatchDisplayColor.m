function imageBatchDisplayColor(D,width,height,nrow,ncol,start,imgNames)
% D: a list of images... each column corresponding to an image
% width/height: the width/height of an image
% nrow/ncol: how many rows/cols you want to have in the layout?
% colorScale: can be '' or 'gray'
    [p,n,junk] = size(D);
    nR = nrow*height;
    nC = ncol*width;
    I = zeros(nR,nC,3);
    rowStart = 1;
    colStart = 1;
    k = 1;
   for i = 1:nrow
       for j = 1:ncol
           indR = (rowStart + (i-1)*height):( i*height);
           indC = (colStart + (j-1)*width):( j*width);
           for q = 1:3
               I(indR,indC,q) = reshape(D(:,k,q),height,width);
           end
           k = k+1;
           if k>n
               break;
           end
       end
       if k>n
           break;
       end
   end
   %if ~exist('r')
   %    r = [min(min(I)),max(max(I))];
   %end
   figure('position',[1404,0,150*ncol,75*nrow]);
   imshow(I);
   
   axis off;        
  
   for j = 1:ncol
       x = [colStart + (j-1)*width,1];
       y = [colStart + (j-1)*width,nR];
       drawLineSeg(x,y);
   end
   for i = 1:nrow
       x = [1,rowStart + (i-1)*height];
       y = [nC,rowStart + (i-1)*height];
       drawLineSeg(x,y);
   end
   
   k = 1;
   for i = 1:nrow
       for j = 1:ncol
           y =  i*height - 8;
           x =  j*width - 8;
           if exist('imgNames')
               
               text(x,y,[num2str(start+k-1),': ',imgNames{k}],'color','red','horizontalAlignment','right');
           else
               text(x,y,num2str(start+k-1),'color','red','horizontalAlignment','center');
           end
           k = k+1;  
           if k > n
               break;
           end
       end
       if k > n
           break;
       end
   end
end
