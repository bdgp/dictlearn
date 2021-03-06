function imageBatchDisplay2(D,width,height,nrow,ncol,imgNames,colorScale,start,showNum,r)
% D: a list of images... each column corresponding to an image
% width/height: the width/height of an image
% nrow/ncol: how many rows/cols you want to have in the layout?
% start: the start of the labeling number
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
   %if ~exist('r')
   %    r = [min(min(I)),max(max(I))];
   %end
   figure('position',[1404,0,150*ncol,75*nrow]);
   if strmatch(colorScale,'gray','exact')
       if ~exist('r')
           r = [min(min(1-I)),max(max(1-I))];
       end
       imagesc(1-I,r);
       colormap(gray);
   else
       if ~exist('r')
           r = [min(min(I)),max(max(I))];
       end
       imagesc(I,r); 
   end
   
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
           y =  i*height - 4;
           x =  j*width - 2;
           if showNum
               text(x,y,[num2str(start+k-1),': ',imgNames{k}],'color','red','horizontalAlignment','right');
           else
               text(x,y,imgNames{k},'color','red','horizontalAlignment','right');
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
