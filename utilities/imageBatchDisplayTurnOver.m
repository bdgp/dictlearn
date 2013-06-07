function imageBatchDisplayTurnOver(D,width,height,nrow,ncol,imgNames)
% D: a list of images
% width/height: the width and height of an image
% nrow/ncol: number of rows/cols of the layout
% imgNames: the names of the images
numImg = size(D,2);
numPerPage = nrow*ncol
numPage = floor(numImg/numPerPage);
for i = 1:(numPage+1)
    indTemp = ((i-1)*numPerPage + 1):min((i*numPerPage),numImg);
    if size(D,3) == 1
        imageBatchDisplay2(D(:,indTemp),width,height,nrow,ncol,imgNames(indTemp),'gray',indTemp(1));
    else
        imageBatchDisplayColor(D(:,indTemp,:),width,height,nrow,ncol,indTemp(1),imgNames(indTemp));
    end 
    
end
