source('loadData.R')
numPatterns = c(6,13,25);
readDict = TRUE;

if (readDict){

   dictPath = paste('../32by16/randomStart/K=',6,'/bestDict.mat',sep='');
   dictFile = readMat(dictPath);
   D = dictFile$Dbest;
   Dstd0 = D;
   for (i in 1:dim(D)[2]){
     Dstd0[,i] = D[,i]/max(D[,i]);
     }
   
   dictPath = paste('../32by16/randomStart/K=',13,'/bestDict.mat',sep='');
   dictFile = readMat(dictPath);
   D = dictFile$Dbest;
   Dstd1 = D;
   for (i in 1:dim(D)[2]){
     Dstd1[,i] = D[,i]/max(D[,i]);
     }

   dictPath = paste('../32by16/randomStart/K=',25,'/bestDict.mat',sep='');
   dictFile = readMat(dictPath);
   D = dictFile$Dbest;
   Dstd2 = D;
   for (i in 1:dim(D)[2]){
     Dstd2[,i] = D[,i]/max(D[,i]);
     }
     
}

lambda1 = 0;
idx0 = reorderDict(Dstd0,template = template);
idx1 = reorderDict(Dstd1,template = template);
idx2 = reorderDict(Dstd2,template = template);

Dstd0 = Dstd0[,idx0];
Dstd1 = Dstd1[,idx1];
Dstd2 = Dstd2[,idx2];

alpha0 = spams.lasso(X,D=Dstd0,lambda1=lambda1,return_reg_path = FALSE,pos =TRUE,ols=TRUE);
alpha1 = spams.lasso(X,D=Dstd1,lambda1=lambda1,return_reg_path = FALSE,pos =TRUE,ols=TRUE);
alpha2 = spams.lasso(X,D=Dstd2,lambda1=lambda1,return_reg_path = FALSE,pos =TRUE,ols=TRUE);

alpha0 = as.matrix(alpha0);
alpha1 = as.matrix(alpha1);
alpha2 = as.matrix(alpha2);
C1 = cor(t(alpha0),t(alpha1));
C2 = cor(t(alpha1),t(alpha2));
CList = list(C1,C2);

linkInfo = linkDict2Hungary(Dstd1,Dstd2,C2,maxLinks = 3,threshold = 0.5,template = template);

savePath = './compareDict/6vs13vs25alphaCorr.png';
D = cbind(Dstd0,Dstd1,Dstd2)
linkDictK(D,numPatterns,CList,maxLink = 3,threshold = 0.5,savePath = savePath,template =template)

#linkDict2(Dstd1,Dstd2,C2,maxLinks = 3,threshold = 0.5,savePath = savePath,template =template)

pause;



C1 = cor(Dstd0,Dstd1);
C2 = cor(Dstd1,Dstd2);
CList = list(C1,C2);

savePath = './compareDict/6vs13vs25spatialCorr.png';

#compareDict(Dstd1,Dstd2,C=C1,width = 32, height =16,template = template,savePath = savePath)
#numPatterns = c(6,13,25);
linkDictK(D,numPatterns,CList,maxLink = 2,threshold = 0.5,savePath = savePath,template =template)
