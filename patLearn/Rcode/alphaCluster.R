source('loadData.R');
numPatterns = c(25,50,75);
Lambda1 = c(1,2,5,10);
readDict = TRUE;
path0 = './dictFit/';
imgIdx = 1:500;

for (k in 1:length(numPatterns)){

if (readDict){
   K = numPatterns[k];
   dictPath = paste('../32by16/randomStart/K=',as.character(K),'/bestDict.mat',sep='');
   dictFile = readMat(dictPath);
   D = dictFile$Dbest;
   Dstd = D;
   for (i in 1:dim(D)[2]){
     Dstd[,i] = D[,i]/max(D[,i]);
     }
}


path1 = paste(path0,'K=',K,'/',sep='');
for (l in 1:length(Lambda1)){
lambda1 = Lambda1[l];
path2 = paste(path1, 'lambda=',lambda1,'/',sep='');
alpha = spams.lasso(X,D=Dstd,lambda1=lambda1,return_reg_path = FALSE,pos =TRUE,ols=TRUE);
alpha = as.matrix(alpha);

numNonZeros = colSums(abs(alpha)>1e-6);
dir.create(path2,recursive = T);



#distImg = dist(t(alpha[,imgIdx])); # Euclidean distance
distImg = as.dist(1 - cor(alpha[,imgIdx],method = 'spearman'));
hr <- hclust(distImg, method="single"); # h-clustering 
savePath = paste(path2,'alphaHeirClusterRankCorrTree.png',sep='');
png(savePath,height=1000,width =1000);
titleText = paste('K=',K,'lambda=',lambda1);
plot(hr,main=titleText);
dev.off();

L0 = hrClusterCoord(hr);

Xtemp = matrix(0,32*16,length(imgIdx));
Xtemp[ind,] = X[,imgIdx];
L0[,2] = 40*L0[,2];
hr$height = 40*hr$height;
out = plotImgOn2dPlane(L0,Xtemp,width=32,height=16,scale1=30,scale2=1,mergeHeight = hr$height);
savePath = paste(path2,'alphaHeirClusterRankCorr.png',sep='');
png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
text(L[,1],L[,2],geneNames[imgIdx],cex = 1,col='red',pos = 1);
drawDendrogram(hr,out$mh,L);
dev.off();

}
}




doHierCluster=0
if (doHierCluster){
#imgIdx = sample(1:dim(X)[2],500);


distImg = dist(t(X[,imgIdx])); # Euclidean distance
distImg = as.dist(1 - cor(X[,imgIdx],method = 'spearman'));
hr <- hclust(distImg, method="single"); # h-clustering 
#pdf('./TFhclustMinDom.pdf',height = 12,width = 14)
#png('./figures/globalHeirCluster.png',height = 1000,width = 1000)
L0 = hrClusterCoord(hr);

Xtemp = matrix(0,32*16,length(imgIdx));
Xtemp[ind,] = X[,imgIdx];

out = plotImgOn2dPlane(L0,Xtemp,width=32,height=16,scale1=30,scale2=1,mergeHeight = hr$height);
savePath = './figures/globalHeirClusterRankCorr.png';
png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
text(L[,1],L[,2],geneNames[imgIdx],cex = 1,col='red',pos = 1);
drawDendrogram(hr,out$mh,L);
dev.off();
}
doKmeans = 0;
if (doKmeans){
km = kmeans(t(X),centers = 49, iter.max=1000,nstart = 10);
imageBatchDisplaySaveMemory(t(km$center),template=template,savePath='./figures/globalClusterKmeansUnstd.png')
png('./figures/globalClusterKmeansUnstdHist.png');hist(km$cluster,breaks=49)
dev.off()
}

#plot(hr,cex = 1,lwd=0.5, xlab = ' ', ylab =' ', 
#     main = 'Global hierarchical clustering of all images via single-linkage')
#dev.off();
doHeatmap = 0;
if (doHeatmap){
library(gplots);
 
mycl <- cutree(hr, k = 49) # end up with many clusters
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(mycl)] # Cuts the tree and creates color vector for clusters.

myheatcol <- rainbow(9);
png('./figures/globalHeirClusterHeatmap.png',height = 1000,width = 1000);
C = cor(Xtemp);
plotOrder = hr$order;
C = C[plotOrder,plotOrder];

heatmap.2(C, Rowv=FALSE, Colv= FALSE, #as.dendrogram(hr),
          dendrogram = "none", 
          col=myheatcol, scale="none", density.info="none", key = FALSE,
          trace="none", ColSideColors=mycolhc,zlim=c(0,1),cexCol = 0.2,cexRow = 0.5)

dev.off();
}