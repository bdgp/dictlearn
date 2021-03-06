plotImgOn2dPlane<-function(l,X,width,height,scale1,scale2,mergeHeight = NULL,topMargin = 0, margin = 10){
#l: N by 2 matrix, (pos) integer-valued
#X: N vectorized images;
#width/height: original width/height of the images 
#scale1: scale of the coord
#scale2: scale of the image

drawBox = 1;
boxVal = max(X)/5;

xMin = min(l[,1]);
xMax = max(l[,1]);
yMin = min(l[,2]);
yMax = max(l[,2]);

L = l; # new coordinate
L[,1] = (l[,1]-xMin)*scale1+height*scale2 + margin;
L[,2] = (l[,2]-yMin)*scale1+width*scale2 + margin;
if (!is.null(mergeHeight)){
   mergeHeight = (mergeHeight-yMin)*scale1+width*scale2 + margin;
}
L = floor(L);

img = matrix(0,nrow = max(L[,2])+margin+topMargin, ncol = max(L[,1])+margin);
imgTemp2 = imresize(matrix(X[,1],nrow = height,ncol = width),scale2);

nrow = dim(imgTemp2)[1];
ncol = dim(imgTemp2)[2];

for (i in 1:dim(L)[1]){
    imgTemp1 = matrix(X[,i],nrow = height,ncol = width);
    imgTemp2 = imresize(imgTemp1,scale2);
     
    if (drawBox){
        imgTemp2[1,] = boxVal;
        imgTemp2[,1] = boxVal;
        imgTemp2[dim(imgTemp2)[1],] = boxVal;
        imgTemp2[,dim(imgTemp2)[2]] = boxVal;
    }
    
    if ((nrow%%2) == 0) {
        rowInd = (L[i,2]- nrow/2+1):(L[i,2]+nrow/2);
    }else{
        rowInd = (L[i,2]- floor(nrow/2)):(L[i,2]+floor(nrow/2));
    }
    # what are you doing here??
    rowInd = rowInd[seq(from = length(rowInd), to = 1, by = -1)];
    if ((ncol%%2) == 0){
        colInd = (L[i,1]- ncol/2+1):(L[i,1]+ncol/2);
    }else{
        colInd = (L[i,1]- floor(ncol/2)):(L[i,1]+floor(ncol/2));
    }
    img[rowInd,colInd] = imgTemp2;
        
   
}
nr = dim(img)[1];
nc = dim(img)[2];

x = L[,1]/nc;
y = L[,2]/nr;

L[,1] = x;
L[,2] = y;

if (!is.null(mergeHeight)){
   mergeHeight = mergeHeight/nr + (height/2)/nr;
   return(list(img = img, L = L, mh = mergeHeight));

}

return(list(img = img, L = L));
}