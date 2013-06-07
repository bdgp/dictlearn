D = Dnew[,1:(dim(Dnew)[2]-2)];
L = patternCentroid(D,template=template);
L = L[1:dim(D)[2],];
K = dim(L)[1];
Dtemp = matrix(0,nrow = 32*16,ncol = dim(D)[2]);
Dtemp[ind,] = D;

out = plotImgOn2dPlane(L,Dtemp,width=32,height=16,scale1=30,scale2=1,margin = 30);
savePath = './patCluster/patOnPlaneTF.png'
png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
#text(L[,1],L[,2],1:dim(D)[2],cex = 1,col='red',pos = 1);

numImgPat = rep(0,dim(L)[1]);
for (i in 1:dim(L)[1]){
    numImgPat[i] = length(unique(geneNames[patImgIdx[[i]]]));
}
cex = 40*numImgPat/max(numImgPat);
points(L[,1],L[,2],pch=16,col='#0000FF50',cex = cex);

colRainbow = rainbow(K*(K-1)/2,alpha = .3);
countCommonGenes = matrix(0,nrow=K,ncol=K)
q = 1;
for (i in 1:(K-1)){
    for (j in (i+1):K){
    	intTemp = intersect(patImgIdx[[i]],patImgIdx[[j]]);
	countCommonGenes[i,j] = length(unique(geneNames[intTemp]));
	if (countCommonGenes[i,j]>5){
	   segments(L[i,1],L[i,2],L[j,1],L[j,2],lwd = (length(intTemp)/8)^2,col=colRainbow[(i-1)*K+j]);
	   q = q+1;}
	}
}
text(L[,1],L[,2],numImgPat,pch=16,col='black',cex = 2);

dev.off();pause


stg = 2;
organ = 'FoGut';
organStageIdx = (which(OSLabel == organ)-1)*5+stg;
geneExpTemp = tfExp[,organStageIdx];
numTfExp = rep(0,dim(L)[1]);
for (i in 1:dim(L)[1]){
    idxTemp = unique(geneNames[patImgIdx[[i]]]);
    sumTemp = 0;
    for (j in 1:length(idxTemp)){
    	tfInd = which(tfNames == idxTemp[j]);
	sumTemp = sumTemp + geneExpTemp[tfInd];
    }
    if (length(sumTemp) == 0){
       numTfExp[i] = 0;
       }else{
       numTfExp[i] = sumTemp;
   }
}
cex = 20*numTfExp/numImgPat;
points(L[,1],L[,2],pch=17,col='#00FF0080',cex = cex);

dev.off()