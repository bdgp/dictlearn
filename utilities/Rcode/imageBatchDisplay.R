imageBatchDisplay<-function(D,width=32,height=16,nrow=7,ncol=7,colorScale = 'gray',imgNames = NULL, savePath = NULL, noNumber = FALSE){
# D: a list of images... each column corresponding to an image
# width/height: the width/height of an image
# nrow/ncol: how many rows/cols you want to have in the layout?
# colorScale: can be '' or 'gray'
   d = dim(D);
   p = d[1]; n = d[2];
   nR = nrow*height;
   nC = ncol*width;
   I = matrix(0,nrow = nR,ncol = nC);
   rowStart = 1;
   colStart = 1;
   k = 1;
   for (i in 1:nrow){
       for (j in 1:ncol){
           indR = (rowStart + (i-1)*height):( i*height);
           indC = (colStart + (j-1)*width):( j*width);
           I[indR,indC] = matrix(D[,k],nrow = height,ncol = width);
           k = k+1;
           if (k>n){
               break;
	   }
        }
  
       if (k>n){
           break;
	   }
   }
   
   imgTemp = t(I);
   indTemp = seq(from = dim(imgTemp)[2], to = 1, by = -1);
   imgTemp = imgTemp[,indTemp];
   if (colorScale == 'gray'){
      colorScale1 = gray(seq(from = 1, to = 0, by = -0.05));
      textCol = 'red'
   }else if (colorScale == 'blueRed'){
   	 colorScale1 = blue2red(50,'ff');
	 colorScale1 = colorScale1[50:1];
	 textCol = 'black';   
   }
   #colorScale = gray(seq(from = 0, to = 1, by = 0.05));
   oldMargin = par('mar');
   if (is.null(savePath)){
      X11(height = nR/30, width = nC/30);
      par(mar = c(0.5,0.5,0.5,0.5));
    }
   else{
      png(file = savePath,height = 3*nR,width = 3*nC);
      par(mar = c(.5,.5,.5,.5));
  }
   if (colorScale == 'gray'){
  
      image(imgTemp,col = colorScale1, xaxt = 'n', yaxt = 'n');
   }else if (colorScale == 'blueRed'){
        lowerV = min(D);
	upperV = max(D);
  	r = max(abs(lowerV),upperV);  
	
      	image(imgTemp,col = colorScale1, xaxt = 'n', yaxt = 'n',zlim =c(-r,r));
   }
   
   #return(I);   
   par(mar = oldMargin);
   
   
   nRow = dim(imgTemp)[1];
   nCol = dim(imgTemp)[2];
   
   horiLines = (seq(from = 1,to = nCol,by = height)-0.5)/nCol;
   vertLines = (seq(from = 1,to = nRow,by = width)-0.5)/nRow;
   
   horiLines[1] = NA;
   vertLines[1] = NA;

   abline(h = horiLines);
   abline(v = vertLines);

    
   k = 1;
   xCoord = NULL;
   yCoord = NULL;
   charTemp = NULL;
   for (i in 1:nrow){
       for (j in 1:ncol){
           y =  1 - (i*height - 3)/nCol;
           x =  (j*width)/nRow;
	   xCoord = c(xCoord,x);
	   yCoord = c(yCoord,y);
   	   if (!is.null(imgNames)){
	      if (noNumber){
	      charTemp = c(charTemp,imgNames[k]);
	      }else{
	      charTemp = c(charTemp, paste(k,': ',imgNames[k],sep = ''));		}
	   }else{
	      if (noNumber){
	      	 charTemp = c(charTemp,' ');
		}else{
	      charTemp = c(charTemp, k);}
	   }
           k = k+1;  
           if (k > n){
               break;
           }
       }
       if (k > n){
           break;
	   }
       
   }
   #charTemp = as.character(1:min(n,nrow*ncol));
   text(xCoord,yCoord,charTemp,col=textCol,cex = 0.8,pos = 2);
   if (!is.null(savePath)){
      dev.off();
     }
}
