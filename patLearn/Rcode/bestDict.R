source('loadData.R');
#Xstd = imageStd(X,'quantile',.5);
numPatterns = seq(from=1,to=100,by=2);
numReplicates = 100;

lambda1 = 0;
gamma1 = 0;
posD = TRUE;
posAlpha =TRUE;
iter = 300;
modeD = 'L2';
mode = 'PENALTY';
batchSize = 1024;

set.seed(215);
path0 = '../32by16/noiseStart/'
for (k in 1:length(numPatterns)){
    K = numPatterns[k];
    path1 = paste(path0,'K=',K,sep='');
    dir.create(path1,recursive = T);

    for (q in 1:numReplicates){
    	
	D0 = matrix(runif(dim(X)[1]*K),dim(X)[1],K);
	D = spams.trainDL(X,K=K,D=D0,lambda1=lambda1,gamma1=gamma1,posD=posD,posAlpha=posAlpha,iter=iter,modeD=modeD,mode=mode,batchsize=batchSize,return_model = FALSE);

	path2 = './rOutput/stdMedian.png';
imageBatchDisplaySaveMemory(D,width=32,height=16,nrow=8,ncol=10,template=template,paintBackground = TRUE,savePath = savePath);
