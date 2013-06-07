source('loadData.R')
numPatterns = c(6,13,25,75);
Lambda1 = c(1);
readDict = TRUE;
path0 = './dictFit/unStd/lassoMean/';


dictPath = paste('../32by16/randomStart/K=1/bestDict.mat',sep='');
dictFile = readMat(dictPath);
D = dictFile$Dbest;
Dmean = D/max(D);

for (k in 1:length(numPatterns)){

if (readDict){
   K = numPatterns[k];
   K = 47
   dictPath = paste('../32by16/randomStart/K=',as.character(K),'/bestDict.mat',sep='');
   dictFile = readMat(dictPath);
   D = dictFile$Dbest;
   Dstd = D;
   for (i in 1:dim(D)[2]){
     Dstd[,i] = D[,i]/max(D[,i]);
     }
   idx0 = reorderDict(Dstd,template = template);
   Dstd = Dstd[,idx0];
}


path1 = paste(path0,'K=',K,'/',sep='');
for (l in 1:length(Lambda1)){
lambda1 = Lambda1[l];
lambda1 = 0;
path2 = paste(path1, 'lambda=',lambda1,'/',sep='');
dir.create(path2,recursive = T);


doLasso = TRUE;
doLassoIntercept = FALSE;
doLassoMean = FALSE;
doLassoDemeanRes = FALSE;
doLassoDemeanResPat = FALSE;

if (doLasso){
   alpha = spams.lasso(X,D=Dstd,lambda1=lambda1, pos = 1,return_reg_path = FALSE,ols=FALSE);
   
}
pause
if (doLassoDemeanResPat){
   Y = X;
   for (ii in 1:ncol(X)){
       Y[,ii] = X[,ii]-mean(X[,ii]);
   }
   for (ii in 1:ncol(Dstd)){
       Dstd[,ii] = Dstd[,ii] - mean(Dstd[,ii]);
   }
   alpha = spams.lasso(Y,D=Dstd,lambda1=lambda1,return_reg_path = FALSE,ols=FALSE);
   
}

if (doLassoDemeanRes){
   Y = X;
   for (ii in 1:ncol(X)){
       Y[,ii] = X[,ii]-mean(X[,ii]);
   }
   alpha = spams.lasso(Y,D=Dstd,lambda1=lambda1,return_reg_path = FALSE,ols=FALSE);
   
}
if (doLassoIntercept){   
   DwithInt = cbind(Dstd,rep(1,nrow(Dstd)));
   W = matrix(1,nrow = ncol(Dstd),ncol = ncol(X));
   W = rbind(W,rep(1e-6,ncol(X)));
   # the intercept is not penalized
   alpha = spams.lassoWeighted(X,D=DwithInt,W = W, lambda1=lambda1);
   alpha = alpha[1:(nrow(alpha)-1),];
}
if (doLassoMean){   
   DwithMean = cbind(Dstd,Dmean);
   W = matrix(1,nrow = ncol(Dstd),ncol = ncol(X));
   W = rbind(W,rep(1e-6,ncol(X)));
   # the intercept is not penalized
   alpha = spams.lassoWeighted(X,D=DwithMean,W = W, lambda1=lambda1);
   Dstd = DwithMean;
   
}
pause

doOLS = TRUE;
if (doOLS){
   if (doLassoDemeanRes|doLassoDemeanResPat){
      result = olsLassoFit(Y,Dstd,as.matrix(alpha));
    }else{
    result = olsLassoFit(X,Dstd,as.matrix(alpha));
    }
   alpha = result$alphaOLS;
   Dnew = result$D;
   numNonZeros = colSums(abs(alpha)>1e-6);
}

png(paste(path2,'nonZeroAlphaHist.png',sep=''));
hist(numNonZeros,breaks = K, main='histogram of number of nonzero coefficients per image', xlab = 'number of nonzeros');
abline(v = mean(numNonZeros),col='red');
dev.off();

png(paste(path2,'rsquareHist.png',sep=''));
hist(result$rSq,breaks = 100, main='histogram of r-squared', xlab = 'r-squared');
abline(v = mean(result$rSq),col='red');
dev.off();

imgIdxManyPat = sort.int(numNonZeros,decreasing = FALSE, index.return = TRUE)$ix;
path3 = paste(path2,'fewPatternsImg/',sep='');
dir.create(path3,recursive = T);

for (i in 1:20){
    imgIdx0 = ((i-1)*10+1):(i*10);
    imgIdx = imgIdxManyPat[imgIdx0];
    
    savePath = paste(path3,imgIdx0[1],'to',imgIdx0[10],'.png',sep = '');
    if (doLassoDemeanRes|doLassoDemeanResPat){ 
       displayFit(Y[,imgIdx],Dnew,alpha[,imgIdx],imgNames=geneNames[imgIdx],template = template,width = 32, height =16,scale =1, maxCovariates = 10, savePath = savePath);
    }else{
       displayFit(X[,imgIdx],Dnew,alpha[,imgIdx],imgNames=geneNames[imgIdx],template = template,width = 32, height =16,scale =1, maxCovariates = 10, savePath = savePath);
       }
}
pause


for (i in 1:20){
    imgIdx = ((i-1)*10+1):(i*10);

    savePath = paste(path2,imgIdx[1],'to',imgIdx[10],'.png',sep = '');
    if (doLassoDemeanRes|doLassoDemeanResPat){ 
       displayFit(Y[,imgIdx],Dnew,alpha[,imgIdx],imgNames=geneNames[imgIdx],template = template,width = 32, height =16,scale =1, maxCovariates = 10, savePath = savePath);
    }else{
       displayFit(X[,imgIdx],Dnew,alpha[,imgIdx],imgNames=geneNames[imgIdx],template = template,width = 32, height =16,scale =1, maxCovariates = 10, savePath = savePath);
       }
}
pause
}
}
    

#imgIdx = 1:10;
#displayFit(X[,imgIdx],Dstd,alpha[,imgIdx],imgNames=geneNames[imgIdx],template = template,width = 32, height =16,scale =1, maxCovariates = 5, savePath = '../rcodeTest/fitTest1.png');


