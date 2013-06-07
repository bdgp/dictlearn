
if (!exists('Load')){
   source('./loadFunctions.R');
   Load = 1;
}

if (!exists('readData')){
   # load the matlab file that contains the gene names for the images
   path0 = "/users/siqi/fruitfly/image/osDict/principalPatternPaper/dictLearn/bestDict/";
   dataFile0 = readMat(paste(path0,'geneNames.mat',sep = ''));
   geneNames = as.vector(unlist(dataFile0$geneName));

   # load the image data 32 by 16 with template.
   dataPath = paste(path0,'32by16/data.mat',sep = '');
   dataFile = readMat(dataPath);
   template = dataFile$template;
   
   X = dataFile$X;
   X[X<0] = 0;
   X[X>1] = 1;
   colnames(X)<-geneNames;
   ind = which(template[,,1]==1);
   
   # Now load the TF annotation data.
   dataFile2 = readMat('~/fruitfly/annotation/TF/annotDataFinal.mat');
   tfNames = as.vector(unlist(dataFile2$gsym));
   tfExp = dataFile2$gmat;
   tfDom = dataFile2$gdom;
   
   tfInd = NULL
   for (i in 1:length(tfNames)){
       indTemp = which(geneNames == tfNames[i]);
       tfInd = c(tfInd,indTemp);
       } 

   numOS = 16;
   OSLabel = c('VisualPr','CNS','Ect/Epi','SalGl','Tracheal','FoGut','SNS',
            'Endo/Midgut','HiGut','Meso/Muscle','Endocrine/Heart',
            'Blood/Fat','Pole/Germ cell','Extraemb','imagPR','PNS');
	    OSLabelStage = rep('NULL',numOS*5);
  for (i in 1:numOS){
    for (j in 1:5){
       labelTemp = paste(OSLabel[i],'_',j,sep = '')
       OSLabelStage[5*(i-1)+j] = labelTemp;
     }
    }
   colnames(tfExp) <- OSLabelStage
    

   readData = 1;

}

stdMax = 0;
if (stdMax){
   Xstd = X;
   for (i in 1:dim(X)[2]){
       Xstd[,i] = X[,i]/max(X[,i]);
       }
}

