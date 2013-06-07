#set.seed(215);
Xstd = X;
K = 25;
lambda1 = 0;
gamma1 = 0;
posD = TRUE;
posAlpha =TRUE;
iter = 500;
modeD = 'L2';#'L1L2';
mode = 'PENALTY';
batchSize = 1024;


imageBatchDisplaySaveMemory(Dstd,width=32,height=16,template=template,paintBackground = TRUE,colorScale ='blueRed',nrow = 7,ncol=7,savePath = './compareDict/testImg/7.png');

alpha = spams.lasso(X,D=Dstd,lambda1=0,pos = TRUE);
alpha = as.matrix(alpha)

alphaD0 = matrix(runif(dim(alpha)[1]*K),dim(alpha)[1],K);
alphaD = spams.trainDL(alpha,K=K,D = alphaD0,lambda1=lambda1,gamma1=gamma1,posD=posD,posAlpha=posAlpha,iter=iter,modeD=modeD,mode=mode,batchsize=batchSize,return_model = FALSE);


imageBatchDisplaySaveMemory(Dstd%*%alphaD,width=32,height=16,template=template,paintBackground = TRUE,colorScale ='blueRed',nrow = 5,ncol=5,savePath = './compareDict/testImg/8.png');

K = 12
alpha = as.matrix(alpha);
alpha2 = spams.lasso(alpha,D=alphaD,lambda1=lambda1,return_reg_path = FALSE,ols=FALSE,pos =TRUE);
alpha2 = as.matrix(alpha2);
alphaD0 = matrix(runif(dim(alpha2)[1]*K),dim(alpha2)[1],K);
alphaDD = spams.trainDL(alpha2,K=K,D = alphaD0,lambda1=lambda1,gamma1=gamma1,posD=posD,posAlpha=posAlpha,iter=iter,modeD=modeD,mode=mode,batchsize=batchSize,return_model = FALSE);
imageBatchDisplaySaveMemory(Dstd%*%alphaD%*%alphaDD,width=32,height=16,template=template,paintBackground = TRUE,colorScale ='blueRed',nrow = 5,ncol=5,savePath = './compareDict/testImg/9.png');


pause

D0 = matrix(runif(dim(X)[1]*K),dim(X)[1],K);
D = spams.trainDL(Xstd,K=K,D=D0,lambda1=lambda1,gamma1=gamma1,posD=posD,posAlpha=posAlpha,iter=iter,modeD=modeD,mode=mode,batchsize=batchSize,return_model = FALSE);

#savePath = './rOutput/stdMedian.png';
Dstd = D;
for (j in 1:dim(D)[2]){
    Dstd[,j] = D[,j]/max(D[,j]);
}
idx0 = reorderDict(Dstd,template = template);
Dstd = Dstd[,idx0];

imageBatchDisplaySaveMemory(Dstd,width=32,height=16,template=template,paintBackground = TRUE,colorScale ='blueRed',nrow = 8,ncol=8);
alpha = spams.lasso(X,D,lambda1=lambda1,return_reg_path = FALSE,pos =FALSE,ols=FALSE);
alpha = as.matrix(alpha);
numNonzero = colSums(abs(alpha)>0.0001);
X11();hist(numNonzero,breaks=K);
abline(v=mean(numNonzero),col='red');


pause

pause;

lambda1 = 2;
Dstd = D;
for (i in 1:dim(D)[2]){
    Dstd[,i] = D[,i]/max(D[,i]);
}

alpha = spams.lasso(X,D=Dstd,lambda1=lambda1,return_reg_path = FALSE,pos =TRUE,ols=TRUE);

imgIdx = 1:10;
displayFit(X[,imgIdx],Dstd,alpha[,imgIdx],imgNames=geneNames[imgIdx],template = template,width = 32, height =16,scale =1, maxCovariates = 5, savePath = '../rcodeTest/fitTest2.png');


#imageBatchDisplaySaveMemory(D,width=32,height=16,nrow=7,ncol=7,template=template,paintBackground = TRUE);
numNonzero = colSums(alpha>0.0001);
X11();hist(numNonzero,breaks=K);
#imageBatchDisplaySaveMemory(D%*%as.matrix(alpha),width=32,height=16,nrow=7,ncol=7,template=template);
