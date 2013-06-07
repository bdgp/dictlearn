numPat = dim(Dnew)[2]-2;#dim(DmultiRes)[2] - 1;
patImgIdx = NULL;
maxNumImgPerPat = 2000;
threshold = 0.05;

doTF = TRUE;
if (doTF){
   dataInd = tfInd;
}else{
   dataInd = 1:length(geneNames);
}

for (i in 1:numPat){
    alphaTemp = alpha[i,dataInd];
    temp = sort.int(abs(alphaTemp),index.return=TRUE,decreasing = TRUE);
    idx = temp$ix;
    val = temp$x;
    n1 = min(which(val<threshold))-1;
    if (n1>0){
       nKeep = min(n1,maxNumImgPerPat);
       idx = idx[1:nKeep];
    }else{
       idx = integer(0);
    }
    patImgIdx = c(patImgIdx,list(dataInd[idx]));
}

printPatImg = FALSE
if (printPatImg){
for (i in 1:numPat){
    savePath = paste('./patCluster/patGeneInteractionTFonly/pat',i,sep='');
    imgTemp = as.matrix(X[,patImgIdx[[i]]]);
    aTemp = alpha[i,patImgIdx[[i]]];
    negIdx = which(aTemp<0);
    imgTemp[,negIdx] = -imgTemp[,negIdx];

    dataTemp = cbind(Dnew[,i],imgTemp);
    patName = paste('Pattern',i);
    gnTemp = rep(' ',length(patImgIdx[[i]]));
    gnTemp2 = geneNames[patImgIdx[[i]]];
    for (j in 1:length(gnTemp)){
    	gnTemp[j] = paste(gnTemp2[j],': ',signif(alpha[i,patImgIdx[[i]]][j],2),sep='');
	}
    #namesTemp = c(patName,geneNames[patImgIdx[[i]]]);
    namesTemp = c(patName, gnTemp);
    imageBatchDisplaySaveMemory(dataTemp,template=template,imgNames = namesTemp, savePath = savePath, nrow = 12,ncol = 10,noNumber = TRUE,colorScale = 'blueRed');

}   
}