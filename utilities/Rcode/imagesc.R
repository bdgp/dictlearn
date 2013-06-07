imagesc <- function(I,colorScale = 'gray',savePath = NULL, invert = TRUE){
imgTemp = t(I);
if (invert){
indTemp = seq(from = dim(imgTemp)[2], to = 1, by = -1);
imgTemp = imgTemp[,indTemp];
}
nR = dim(I)[1];
nC = dim(I)[2];


if (colorScale == 'gray'){
   colorScale = gray(seq(from = 1, to = 0, by = -0.05));
   }else if (colorScale == 'heatmap'){
   colorScale = rainbow(100);   
}
oldMargin = par('mar');
if (is.null(savePath)){
    #X11(height = nR/100, width = nC/100);
    par(mar = c(0.5,0.5,0.5,0.5));
}else{
    png(file = savePath,height = nR,width = nC);
    par(mar = c(.5,.5,.5,.5));
}

image(imgTemp,col = colorScale, xaxt = 'n', yaxt = 'n');    
par(mar = oldMargin);

if (!is.null(savePath)){
    dev.off();
}


  
} 