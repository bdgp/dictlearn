imageBatchDisplaySaveMemory<-function(D,width=32,height=16,nrow=7,ncol=7,colorScale = 'gray',template,paintBackground = FALSE,imgNames = NULL, savePath = NULL,noNumber = FALSE){
# D: a list of images... each column corresponding to an image
# width/height: the width/height of an image
# nrow/ncol: how many rows/cols you want to have in the layout?
# colorScale: can be '' or 'gray'
ind = which(template[,,1]==1);
Dfull = matrix(0,nrow = height*width,ncol = dim(D)[2]);
Dfull[ind,] = D;
if (paintBackground){
      #r = quantile(D,0.94);	
      Dfull[-ind,] = max(D)/3;
}

I = imageBatchDisplay(Dfull, width, height, nrow,ncol,colorScale,imgNames,savePath,noNumber = noNumber); 
return(I);
}