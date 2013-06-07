# an analysis of standardization schemes
# max, L1, L2, z-stat+;
width = 32; height = 16;
X[X<0]=0;
X[X>1]=1;

XstdQuarter = XstdMedian = XstdMean = XstdL2 = XstdL1 = XstdMax = X;
for (i in 1:dim(X)[2]){
    XstdMax[,i] = X[,i]/max(X[,i]);
    XstdL1[,i] = X[,i]/sum(X[,i]);
    XstdL2[,i] = X[,i]/sqrt(sum(X[,i]^2));

    temp = (X[,i] - mean(X[,i])); temp[temp<0]=0;
    XstdMean[,i] = temp;
    
    temp = (X[,i] - median(X[,i])); temp[temp<0]=0;    
    XstdMedian[,i] = temp;
    
    temp = (X[,i] - quantile(X[,i],0.25)); temp[temp<0]=0;    
    XstdQuarter[,i] = temp;


}    


for (i in 1:160){
    imgIdx = ((i-1)*10+1):(i*10);

#imgIdx = 1:20;
Xtemp = X[,imgIdx];    
Y = matrix(1/3,width*height,(1+4)*length(imgIdx));

k = 1;
imgNames = rep(' ',length(imgIdx));
for (i in 1:length(imgIdx)){
    imgNames[k] = geneNames[imgIdx[i]];
    Y[ind,k] = X[,imgIdx[i]]; k = k+1;
    Y[ind,k] = XstdMax[,imgIdx[i]]; k = k+1;
    Y[ind,k] = XstdMean[,imgIdx[i]]; k = k+1;
    Y[ind,k] = XstdMedian[,imgIdx[i]]; k = k+1;
    Y[ind,k] = XstdQuarter[,imgIdx[i]]; k = k+1;
}     

savePath = paste('./std/',imgIdx[1],'to',imgIdx[10],'.png',sep = '');
imageBatchDisplay(Y,width=32,height=16,nrow=length(imgIdx),ncol=(1+4),colorScale = 'gray',imgNames = imgNames, savePath = savePath,noNumber= TRUE);
}
 