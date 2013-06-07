
readData = TRUE
if (readData){
K = 47;
dataFile0 = readMat('../dataTemp.mat');
geneNames = as.vector(unlist(dataFile0$geneName));

dictPath = paste('./32by16/randomStart/K=',as.character(K),'/bestDict.mat',sep='');
dataPath = paste('./32by16/data.mat');
dataFile = readMat(dataPath);
dictFile = readMat(dictPath);
template = dataFile$template;
D = dictFile$Dbest;
X = dataFile$X;
}
ind = which(template[,,1]==1);
imgIdx = 1:50;
Xtemp = matrix(0,32*16,length(imgIdx));
Xtemp[ind,] = X[,imgIdx];
#l = matrix(0,nrow = length(imgIdx),ncol = 2);
#l[,1] = rnorm(length(imgIdx));
#l[,2] = rnorm(length(imgIdx));
out = plotImgOn2dPlane(l,Xtemp,width=32,height=16,scale1=100,scale2=1);

savePath = './rcodeTest/layout.png';
png(file = savePath,height = dim(out$img)[2],width = dim(out$img)[1]);
imagesc(out$img);
L = out$L;
text(L[,1]+18/dim(out$img)[1],L[,2]-4/dim(out$img)[2],geneNames[imgIdx],cex = 0.5,col='red',pos = 2);
dev.off();

#imageBatchDisplaySaveMemory(X[,imgIdx],width=32,height=16,nrow=7,ncol=7,template=template,imgNames = geneNames[imgIdx],savePath ='./rcodeTest/temp.png');
#imageBatchDisplaySaveMemory(D,width=32,height=16,nrow=7,ncol=7,template=template,paintBackground = TRUE);
