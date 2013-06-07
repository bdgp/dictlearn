numPatterns = c(6,13,25,47,1);
readDict = TRUE;
if (readDict){
   DmultiRes = matrix(0,nrow= length(ind),ncol=sum(numPatterns));
   whichDict = rep(0,sum(numPatterns));
   patternNames = rep(' ',sum(numPatterns)+1);
   q = 1;
   for (k in 1:length(numPatterns)){
       K = numPatterns[k];
       dictPath = paste('../32by16/randomStart/K=',K,'/bestDict.mat',sep='');
       dictFile = readMat(dictPath);
       D = dictFile$Dbest;
       Dstd = D;
       for (i in 1:dim(D)[2]){
           Dstd[,i] = D[,i]/max(D[,i]);
	   patternNames[q+i-1] = paste('K',K,'-p',i,sep='');
       }
       idx0 = reorderDict(Dstd,template = template);
       Dstd = Dstd[,idx0];      
       
       idx = q:(q+K-1);
       DmultiRes[,idx] = Dstd;
       whichDict[idx] = k;
       q = q + K;
   }
   patternNames[q] = 'b_0';
}


lambda1 = 1.2;
pruneDict = FALSE;

n0 = length(patternNames);       
if (pruneDict){
   keepDictInd = c(keepDictIdx,n0-1,n0);
   Dnow = DmultiRes[,c(keepDictIdx,n0-1)]; # not include the intercept
   alpha = spams.lasso(X,D=Dnow,lambda1=lambda1,return_reg_path = FALSE,ols=FALSE);
   
}else{
   keepDictInd = 1:n0;
   alpha = spams.lasso(X,D=DmultiRes,lambda1=lambda1,return_reg_path = FALSE,ols=FALSE);
   result = olsLassoFit(X,DmultiRes,as.matrix(alpha));
   Dnow = DmultiRes;
}

doOLS = TRUE;
if (doOLS){
   intercept = TRUE
   result = olsLassoFit(X,Dnow,as.matrix(alpha),intercept = intercept);
   n0 = dim(alpha)[1]; 
   getRidMatTemp = matrix(result$alphaOLS[1:n0,]*alpha<= (-1e-6),nrow=dim(alpha)[1],ncol=dim(alpha)[2]);   
   alpha[getRidMatTemp] = 0; 
   result = olsLassoFit(X,Dnow,as.matrix(alpha),intercept = intercept);
   Dnew = result$D;
   alpha = result$alphaOLS;
}else{
   Dnew = Dnow;
}

complexityAna = FALSE
if (complexityAna){
   colTemp = colSums(abs(alpha));
   #alphaTemp = rep(0,dim(alpha)[2]);
   alphaTemp = matrix(0,nrow = length(numPatterns),ncol = dim(alpha)[2]);
   for (i in 1:dim(alpha)[2]){
       for (k in 1:length(numPatterns)){
              dictIdx = which(whichDict == k);
      	      alphaTemp[k,i] = sum(abs(alpha[dictIdx,i]))
	    }
   }
   #sortIndTemp = (sort.int(alphaTemp,decreasing=TRUE,index.return=TRUE))$ix;
   #imageBatchDisplaySaveMemory(X[,sortIndTemp],template=template,nrow=9,ncol=9);
}

fitData = FALSE
if (fitData){
      
path0 = './multiResFit/testImg/negPosElim/';
if (!doOLS){
   path0 = './multiResFit/testImg/nonOLS/';
   cmd = paste('mkdir ',path0,sep='');
   system(cmd);
}

numNonZeros = colSums(abs(alpha)>1e-6);
png(paste(path0,'nonZeroAlphaHist.png',sep=''));
hist(numNonZeros,breaks = sum(numPatterns), main='histogram of number of nonzero coefficients per image', xlab = 'number of nonzeros');
abline(v = mean(numNonZeros),col='red');
dev.off();

if (doOLS){
png(paste(path0,'rsquareHist.png',sep=''));
hist(result$rSq,breaks = 100, main='histogram of r-squared', xlab = 'r-squared');
abline(v = mean(result$rSq),col='red');
dev.off();
}

for (i in 1:20){
    imgIdx = ((i-1)*10+1):(i*10);

    savePath = paste(path0,imgIdx[1],'to',imgIdx[10],'.png',sep = '');
    displayFit(X[,imgIdx],Dnew,alpha[,imgIdx],imgNames=geneNames[imgIdx],dictNames = patternNames[keepDictInd], template = template,width = 32, height =16,scale =1, maxCovariates = 10, savePath = savePath)

}
}

treeAna = FALSE
if (treeAna){
   geneNamesTemp = c('tll','Doc2','pdm2','pni','CG6464');
   imgIdx = NULL;
   for (i in 1:length(geneNames)){
       idxTemp = which(geneNames == geneNamesTemp[i]);
       imgIdx = c(imgIdx,idxTemp);
   }
   savePath = paste(path0,'tree/exampleImg.png',sep='');
   imageBatchDisplaySaveMemory(X[,imgIdx],template=template,imgNames = geneNames[imgIdx],savePath = savePath);
   for (i in 1:length(imgIdx)){
      coefVal = alpha[1:sum(numPatterns),imgIdx[i]];
      savePath = paste(path0,'tree/',i,geneNames[imgIdx[i]],'.png',sep='');
      linkDictK(DmultiRes,numPatterns,CList,maxLink = 3,threshold = 0.5,savePath = savePath,template =template,coefVal = coefVal);
   }
}

dictUsage = FALSE

if(dictUsage){
    usage = rep(0,sum(numPatterns));
    for (i in 1:length(usage)){
    	usage[i] = sum(abs(alpha[i,])>1e-06);
	}
	pause
	usage = sqrt(10*usage/max(usage));
	savePath = paste(path0,'dictUsageAlphaCountTree.png',sep ='');
      linkDictK(DmultiRes,numPatterns,CList,maxLink = 3,threshold = 0.5,savePath = savePath,template =template,coefVal = usage);


}

