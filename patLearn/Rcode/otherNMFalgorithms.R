methods = list('brunet','lee','nsNMF','offset');
#res = nmf(X,rank=49,method=methods);
for (i in 1:length(methods)){
    D = basis(res[[i]]);
    Dstd = D;
    for (j in 1:dim(D)[2]){
     	Dstd[,j] = D[,j]/max(D[,j]);
    }
    idx0 = reorderDict(Dstd,template = template);
    Dstd = Dstd[,idx0];
    alpha = coef(res[[i]]);
    savePath = paste('./otherNMFalgorithmResults/',methods[i],'49.png',sep='');
    imageBatchDisplaySaveMemory(Dstd,width=32,height=16,template=template,paintBackground = TRUE,colorScale ='blueRed',savePath = savePath);
    png(paste('./otherNMFalgorithmResults/',methods[i],'49alpha.png',sep=''));
    hist(alpha,breaks=100);
    dev.off()
}
    