source('loadData.R')
numPatterns = c(6,13,25,47);
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
}

imageBatchDisplaySaveMemory(DmultiRes,nrow = 10,ncol=10,savePath = './compareDict/multiResDict.png',template =template,imgNames = patternNames)

alphaCorr = FALSE;
if (alphaCorr){
lambda1= 0
alphaList = NULL;
q = 1;
for (i in 1:length(numPatterns)){
    D = DmultiRes[,q:(q+numPatterns[i]-1)];
    alpha = spams.lasso(X,D=D,lambda1=lambda1,return_reg_path = FALSE,pos =TRUE,ols=TRUE);
    alphaList = c(alphaList, list(as.matrix(alpha)));
    q = q+numPatterns[i];
}

CList = NULL;
for (i in 1:(length(numPatterns)-1)){
    alpha1 = alphaList[[i]];
    alpha2 = alphaList[[i+1]];
    C = cor(t(alpha1),t(alpha2));    
    CList = c(CList,list(C));
}
savePath = './compareDict/6vs13vs25vs47alphaCorr.png';

}


spatialCorr = FALSE;
# some linkage patterns are missing due to negligible correlation.
if(spatialCorr){
CList = NULL;
q = 1;
for (i in 1:(length(numPatterns)-1)){
    D1 = DmultiRes[,q:(q+numPatterns[i]-1)];
    d = q+numPatterns[i];
    D2 = DmultiRes[,d:(d+numPatterns[i+1]-1)];
    q = q+numPatterns[i];
    C = cor(D1,D2);    
    CList = c(CList,list(C));
}
savePath = './compareDict/6vs13vs25vs47spatialCorr.png';
}

prunePat = TRUE
if (prunePat){
   threshold = .8
   Dnow = DmultiRes[,1:numPatterns[1]];
   keepDictIdx = 1:numPatterns[1];
   q = 1
   for (i in 1:(length(numPatterns)-1)){       
       d = q+numPatterns[i];
       dictIdx = d:(d+numPatterns[i+1]-1);
       D2 = DmultiRes[,dictIdx];
       q = q+numPatterns[i];
       C = cor(Dnow,D2);
       idx2remove = dictIdx[which(colSums(C >= threshold)>=1)];
       idx2retain = setdiff(dictIdx,idx2remove);
       keepDictIdx = c(keepDictIdx,idx2retain);
       Dnow = cbind(Dnow,DmultiRes[,idx2retain]);
       }
}
   
pause

#compareDict(Dstd1,Dstd2,C=C1,width = 32, height =16,template = template,savePath = savePath)
#linkDict2(Dstd1,Dstd2,C,maxLink = 3,threshold = 0.5,savePath = savePath,template =template)

coefVal = rnorm(sum(numPatterns));
linkDictK(DmultiRes,numPatterns,CList,maxLink = 3,threshold = 0.4,savePath = savePath,template =template,coefVal = NULL);
