compareDict<-function(D1,D2,C,width = 32, height =16,template = template, noValue = FALSE,noNum = FALSE,savePath = NULL){
# D1,D2: the two dictionaries for comparison, D1 will be on the y-axis.
# C: a matrix to be plotted, of size K1 by K2, where K1,K2 are the number of columns of D1 and D2 respectively. below zero values are truncated to zeros.

ind = which(template[,,1]==1);

K1 = dim(D1)[2];
K2 = dim(D2)[2];

Y = matrix(1/3,nrow = width*height,ncol=(K1+1)*(K2+1));
imgNames = rep(' ', (K1+1)*(K2+1));
# fill in the first row with D2
k = 2;
for (i in 2:(K2+1)){
    Y[ind,i] = D2[,i-1];
    imgNames[k] = paste('p',i-1,sep='');
    k = k + 1;
}
for (j in 2:(K1+1)){
    Y[ind,k] = D1[,j-1];
    imgNames[k] = paste('p',j-1,sep='');
    k = k + 1;
    for (i in 2:(K2+1)){
    	Y[,k] = max(C[j-1,i-1],0);
	imgNames[k] = signif(C[j-1,i-1],2);
   	k = k + 1;
    }
}  
    
imageBatchDisplay(Y,width,height,nrow=K1+1,ncol=K2+1,colorScale = 'gray',imgNames = imgNames, savePath = savePath,noNumber= TRUE);

}

linkDictCoord<-function(K1,maxLinks=3,gap = 0.3){
xCoord = yCoord = rep(0, K1 + K1*maxLinks);
xCoord[1:K1] = 1;
xCoord[(K1+1):(K1+K1*maxLinks)] = 2;
yCoord[1:K1] = K1:1;
indTemp = (K1+1):(K1+maxLinks);
for (i in 1:K1){
    
    midpt = yCoord[i];
    if (maxLinks%%2==0){
       startPt = midpt - gap/2 - (maxLinks/2-1)*gap;
    }else{
       startPt = midpt - ((maxLinks-1)/2)*gap;
    }
    addLength = gap*(1:maxLinks-1);
    yCoord[indTemp] = startPt + addLength;
    indTemp = indTemp + maxLinks;
}
L = cbind(xCoord,yCoord);
return(L)

}

linkDictCoordSingle<-function(L0,numLinks=3,gap = 0.3){

xCoord = yCoord = rep(0, numLinks);
xCoord[1:numLinks] = L0[1]+1;

    
midpt = L0[2];
if (numLinks%%2==0){
       startPt = midpt - gap/2 - (numLinks/2-1)*gap;
}else{
       startPt = midpt - ((numLinks-1)/2)*gap;
}
addLength = gap*(1:numLinks-1);
yCoord = startPt + addLength;
yCoord = yCoord[length(yCoord):1];
L = cbind(xCoord,yCoord);
return(L)

}



linkDict2<-function(D1,D2,C,maxLinks = 3,threshold = 0.5,width = 32, height =16,template = template, noValue = FALSE,noNum = FALSE,savePath = NULL){

ind = which(template[,,1]==1);

K1 = dim(D1)[2];
K2 = dim(D2)[2];

Y = matrix(1/3,nrow = width*height,ncol=K1+K1*maxLinks);
imgNames = rep(' ', K1+K1*maxLinks);
L0 = matrix(0, nrow= K1,ncol = 2);
L0[,1] = 1;
L0[,2] = K1:1;

L1 = matrix(0, nrow= K1*maxLinks,ncol = 2);
# fill in the first column with D1
k = 1;
for (i in 1:K1){
    Y[ind,i] = D1[,i];
    imgNames[k] = paste('p',i,sep='');
    k = k + 1;
}
startK = k;

arrowFrom = rep(0,K1*maxLinks);
arrowTo = rep(0,K1*maxLinks);
corrTemp = arrowTo;



Ctemp = matrix(0,nrow=dim(C)[1],ncol=dim(C)[2]);
p = 1;
for (i in 1:K1){
    cTemp = C[i,];
    temp = sort.int(cTemp,decreasing =TRUE, index.return = TRUE);
    indTemp = temp$ix;
    valTemp = temp$x;
    maxIndex = min(which(valTemp<=threshold))-1;
    numLinks = min(maxIndex,maxLinks); 
    Ctemp[i,indTemp[1:numLinks]] = valTemp[1:numLinks];
    ixTemp = p:(p+numLinks-1); 
    
    L1[ixTemp,] = linkDictCoordSingle(L0[i,],numLinks);
    
    p = p+numLinks;

}

L1 = L1[1:(p-1),];
L = rbind(L0,L1);


p = 1;

for (i in 1:K1){
    indNonZeros = which(Ctemp[i,]!=0);
    numLinks = length(indNonZeros);
    Y[ind,startK-1+(1:numLinks)] = D2[,indNonZeros];
    arrowFrom[p-1+(1:numLinks)] = i;
    arrowTo[p-1+(1:numLinks)] = startK-1+(1:numLinks);
    corrTemp[p-1+(1:numLinks)] = Ctemp[i,indNonZeros];
    
    for (j in 1:numLinks){
    	imgNames[startK+j-1] = paste('p',indNonZeros[j],sep='');
	}
    startK = startK + numLinks;
    p = p+numLinks;

}

L[,1] = L[,1]*3;
L[,2] = L[,2]*2;




Y = Y[,1:(startK-1)];
imgNames = imgNames[1:(startK-1)];

out = plotImgOn2dPlane(L,Y,width=32,height=16,scale1=30,scale2=1,margin = 30);

png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
text(L[,1],L[,2],imgNames,cex = 1,col='red',pos = 1);

for (i in 1:length(arrowFrom)){
    x0 = L[arrowFrom[i],1]+ width/(2*dim(out$img)[2]);
    y0 = L[arrowFrom[i],2];
    x1 = L[arrowTo[i],1] - width/(2*dim(out$img)[2]);
    y1 = L[arrowTo[i],2];
    if (arrowTo[i]*arrowFrom[i]==0){break;}
    arrows(x0,y0,x1,y1);
    text((x0+x1)/2,(y0+y1)/2,signif(corrTemp[i],2),col='red');
    
}

dev.off();

}

linkDict2<-function(D1,D2,C,maxLinks = maxLinks,threshold = 0.5,width = 32, height =16,template = template, noValue = FALSE,noNum = FALSE,savePath = NULL){

linkInfo = linkDict2Hungary(D1,D2,C,maxLinks = maxLinks,threshold = threshold, template = template);
ind = which(template[,,1]==1);

K1 = dim(D1)[2];
K2 = dim(D2)[2];

Y = matrix(1/3,nrow = width*height,ncol=K1+K1*maxLinks);
imgNames = rep(' ', K1+K1*maxLinks);
L0 = matrix(0, nrow= K1,ncol = 2);
L0[,1] = 1;
L0[,2] = K1:1;

L1 = matrix(0, nrow= dim(linkInfo)[1],ncol = 2);
# fill in the first column with D1
k = 1;
for (i in 1:K1){
    Y[ind,i] = D1[,i];
    imgNames[k] = paste('p',i,sep='');
    k = k + 1;
}
startK = k;

arrowFrom = rep(0,K1*maxLinks);
arrowTo = rep(0,K1*maxLinks);
corrTemp = arrowTo;



Ctemp = matrix(0,nrow=dim(C)[1],ncol=dim(C)[2]);
p = 1;
for (i in 1:K1){
    indTemp = which(linkInfo[,1] == i);
    numLinks = length(indTemp);
    ixTemp = p:(p+numLinks-1)      
    L1[ixTemp,] = linkDictCoordSingle(L0[i,],numLinks);    
    p = p+numLinks;
}

L1 = L1[1:(p-1),];
L = rbind(L0,L1);


p = 1;

for (i in 1:K1){
    indTemp = which(linkInfo[,1]==i);
    indNonZeros = linkInfo[indTemp,2];
    numLinks = length(indNonZeros);
    Y[ind,startK-1+(1:numLinks)] = D2[,indNonZeros];
    arrowFrom[p-1+(1:numLinks)] = i;
    arrowTo[p-1+(1:numLinks)] = startK-1+(1:numLinks);
    corrTemp[p-1+(1:numLinks)] = C[i,indNonZeros];
    
    for (j in 1:numLinks){
    	imgNames[startK+j-1] = paste('p',indNonZeros[j],sep='');
	}
    startK = startK + numLinks;
    p = p+numLinks;

}

L[,1] = L[,1]*3;
L[,2] = L[,2]*2;




Y = Y[,1:(startK-1)];
imgNames = imgNames[1:(startK-1)];

out = plotImgOn2dPlane(L,Y,width=32,height=16,scale1=30,scale2=1,margin = 30);

png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
text(L[,1],L[,2],imgNames,cex = 1,col='red',pos = 1);

for (i in 1:length(arrowFrom)){
    x0 = L[arrowFrom[i],1]+ width/(2*dim(out$img)[2]);
    y0 = L[arrowFrom[i],2];
    x1 = L[arrowTo[i],1] - width/(2*dim(out$img)[2]);
    y1 = L[arrowTo[i],2];
    if (arrowTo[i]*arrowFrom[i]==0){break;}
    arrows(x0,y0,x1,y1);
    text((x0+x1)/2,(y0+y1)/2,signif(corrTemp[i],2),col='red');
    
}

dev.off();

}

linkDict2Hungary<-function(D1,D2,C,maxLinks = 3,threshold = 0.5,template = template){
# This function applies the hungarian algorithm multiple times to match D1 and D2
# return the matched indices.

K1 = dim(D1)[2];
K2 = dim(D2)[2];

Cost = matrix(0,nrow = K1,ncol = K1*maxLinks);
Cost[,1:K2] = C + 1;

linkInfo = matrix(0,nrow = K1*maxLinks, ncol = 2)
p = 1
unusedInd = 1:(K1*maxLinks);
for (r in 1:maxLinks){
    
    indtemp = as.vector(solve_LSAP(Cost[,unusedInd],maximum = TRUE));
    ind1 = unusedInd[indtemp]; # original index coordinate
    corTemp = rep(0,length(ind1));
   
   
    for (i in 1:length(ind1)){  
    	
    	corTemp[i] = Cost[i,ind1[i]]-1;
	}
	
    keepIdx = which(corTemp>=threshold);
    
    n = length(keepIdx);
    if (n>0){  
    linkInfo[p:(p+n-1),1] = (1:K1)[keepIdx];    
    linkInfo[p:(p+n-1),2] = ind1[keepIdx];
    unusedInd = setdiff(unusedInd, ind1[keepIdx]);
    }
    p = p + n 
}      
linkInfo = linkInfo[1:(p-1),];
indTemp = sort.int(linkInfo[,1],index.return = TRUE)$ix;
linkInfo = linkInfo[indTemp,];

for (i in 1:K1){
    ind1 = which(linkInfo[,1] == i);
    K2Ind = linkInfo[ind1,2];
    ixTemp2 = reorderDict(D2[,K2Ind],template=template);
    linkInfo[ind1,2] = K2Ind[ixTemp2];
}

return(linkInfo);
}


linkDictK<-function(D,numPatterns,CList,maxLinks = 3,threshold = 0.5, width = 32, height =16,template = template, noValue = FALSE,noNum = FALSE,savePath = NULL){

ind = which(template[,,1]==1);
numDict = length(numPatterns);
maxNumImg = sum(numPatterns[1:numDict]*maxLinks) + numPatterns[1];


Y = matrix(1/3,nrow = width*height,ncol=maxNumImg);
imgNames = rep(' ', maxNumImg);

D1 = D[,1:numPatterns[1]];
K1 = dim(D1)[2];

k = 1;
for (i in 1:K1){
     Y[ind,i] = D1[,i];
     imgNames[k] = paste('p',i,sep='');
     k = k + 1;
}
startK = k;

L = matrix(0, nrow= K1,ncol = 2);
L[,1] = 1;
L[,2] = K1:1;


P = 1;
arrowFromList = NULL;
arrowToList = NULL;
corrTempList = NULL;
patternIndList = list(1:K1);
numRows =rep(0,numDict);
numRows[1] = K1;

gaps = c(0.4,0.2,0.1);#0.6*0.6^dictIdx
for (dictIdx in 1:(numDict-1)){
    if (dictIdx == 1){
       indTemp = 1:K1;
    }else{
	indTemp = sum(numPatterns[1:(dictIdx-1)])+(1:numPatterns[dictIdx]);
    }
    D1 = D[,indTemp];
    D2 = D[,sum(numPatterns[1:dictIdx])+(1:numPatterns[dictIdx+1])];
    K1 = dim(D1)[2];
    K2 = dim(D2)[2];


    patternInd = patternIndList[[dictIdx]];

    L0 = L[(startK-numRows[dictIdx]):(startK-1),];
    
    L1 = matrix(0, nrow= K1*maxLinks,ncol = 2);

    arrowFrom = rep(0,K1*maxLinks);
    arrowTo = rep(0,K1*maxLinks);
    corrTemp = arrowTo;
    patternIndNew = arrowTo; 
    C = CList[[dictIdx]];
    Ctemp = matrix(0,nrow=dim(C)[1],ncol=dim(C)[2]);
    p = 1;

    # calculate the coordinates of the second dictionary
    for (j in 1:numRows[dictIdx]){
    	i = patternInd[j];
    	cTemp = C[i,];
    	temp = sort.int(cTemp,decreasing =TRUE, index.return = TRUE);
    	indTemp = temp$ix;
    	valTemp = temp$x;
    	maxIndex = min(which(valTemp<=threshold))-1;
    	numLinks = min(maxIndex,maxLinks); 
    	Ctemp[i,indTemp[1:numLinks]] = valTemp[1:numLinks];
    	ixTemp = p:(p+numLinks-1);     
	L1[ixTemp,] = linkDictCoordSingle(L0[j,],numLinks,gap = gaps[dictIdx]);
	#patternIndNew[ixTemp] = indTemp[1:numLinks];    
	indNZ = which(Ctemp[i,]!=0);
	#print(indNZ)
	ixTemp2 = reorderDict(D2[,indNZ],template=template);
	
	patternIndNew[ixTemp] = indNZ[ixTemp2];
	p = p+numLinks;
    }
    
    patternIndNew = patternIndNew[1:(p-1)];
    numRows[dictIdx+1] = p-1;
    L1 = L1[1:(p-1),];
    L = rbind(L,L1);
    patternIndList = c(patternIndList, list(patternIndNew));
    p = 1;


    # copy the images from the second dictionary
    for (I in 1:numRows[dictIdx]){
    	i = patternInd[I];
	
    	indNonZeros = which(Ctemp[i,]!=0);
	
    	numLinks = length(indNonZeros);
	
	cTemp = Ctemp[i,];
    	temp = sort.int(cTemp,decreasing =TRUE, index.return = TRUE);
    	indTemp = temp$ix;    	
	
	ixTemp2 = reorderDict(D2[,indNonZeros],template=template);
	#patternIndNew[ixTemp] = indNZ[ixTemp2];
    	
    	#Y[ind,startK-1+(1:numLinks)] = D2[,indTemp[1:numLinks]];
	Y[ind,startK-1+(1:numLinks)] = D2[,indNonZeros[ixTemp2]];
	#cat(i,'--->',indNonZeros[ixTemp2],'\n');
	if (dictIdx == 2){
	   #imageBatchDisplay(Y,nrow = 7,ncol = 7);
	   #imageBatchDisplaySaveMemory(D1,template=template);

	   #pause
	   }
	

	if (dictIdx == 1){
	   sumTemp = I;
	   }else{
	   sumTemp = I + sum(numRows[1:(dictIdx-1)]);
	   }
    	arrowFrom[p-1+(1:numLinks)] = sumTemp;
    	arrowTo[p-1+(1:numLinks)] = startK-1+(1:numLinks);
    	corrTemp[p-1+(1:numLinks)] = Ctemp[i,indNonZeros];
    
	for (j in 1:numLinks){
    	    imgNames[startK+j-1] = paste('p',indNonZeros[j],sep='');
	}
    	startK = startK + numLinks;
    	p = p+numLinks;
    }
    arrowFromList = c(arrowFromList,list(arrowFrom));
    arrowToList = c(arrowToList,list(arrowTo));
    corrTempList = c(corrTempList, list(corrTemp));
}

L[,1] = L[,1]*3;
L[,2] = L[,2]*6;
Y = Y[,1:(startK-1)];
imgNames = imgNames[1:(startK-1)];

out = plotImgOn2dPlane(L,Y,width=32,height=16,scale1=30,scale2=.8,margin = 30);

png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
text(L[,1],L[,2],imgNames,cex = 1,col='red',pos = 1);

for (dictIdx in 1:(numDict-1)){
    arrowFrom = arrowFromList[[dictIdx]];
    arrowTo = arrowToList[[dictIdx]];
    corrTemp = corrTempList[[dictIdx]];
    
    for (i in 1:length(arrowFrom)){
    	x0 = L[arrowFrom[i],1]+ width/(2*dim(out$img)[2]);
    	y0 = L[arrowFrom[i],2];
    	x1 = L[arrowTo[i],1] - width/(2*dim(out$img)[2]);
    	y1 = L[arrowTo[i],2];
    	if (arrowTo[i]*arrowFrom[i]==0){break;}
    	   arrows(x0,y0,x1,y1);
    	   text((x0+x1)/2,(y0+y1)/2,signif(corrTemp[i],2),col='red');
    }
}

dev.off();

}



linkDictK<-function(D,numPatterns,CList,maxLinks = 3,threshold = 0.5, width = 32, height =16,template = NULL, noValue = FALSE,noNum = FALSE,savePath = NULL,coefVal = NULL){

if(is.null(template)){
   ind = 1:(width*height);
}

ind = which(template[,,1]==1);
numDict = length(numPatterns);
maxNumImg = sum(numPatterns[1:numDict]*maxLinks) + numPatterns[1];


Y = matrix(1/3,nrow = width*height,ncol=maxNumImg);
imgNames = rep(' ', maxNumImg);

D1 = D[,1:numPatterns[1]];
K1 = dim(D1)[2];

k = 1;
for (i in 1:K1){
     Y[ind,i] = D1[,i];
     imgNames[k] = paste('p',i,sep='');
     k = k + 1;
}
startK = k;

L = matrix(0, nrow= K1,ncol = 2);
L[,1] = 1;
L[,2] = (K1:1)*1.1;


P = 1;
arrowFromList = NULL;
arrowToList = NULL;
corrTempList = NULL;
patternIndList = list(1:K1);
numRows =rep(0,numDict);
numRows[1] = K1;

gaps = c(0.4,0.18,0.08);
for (dictIdx in 1:(numDict-1)){
    if (dictIdx == 1){
       indTemp = 1:K1;
    }else{
	indTemp = sum(numPatterns[1:(dictIdx-1)])+(1:numPatterns[dictIdx]);
    }
    D1 = D[,indTemp];
    D2 = D[,sum(numPatterns[1:dictIdx])+(1:numPatterns[dictIdx+1])];
    K1 = dim(D1)[2];
    K2 = dim(D2)[2];

    
    patternInd = patternIndList[[dictIdx]];

    L0 = L[(startK-numRows[dictIdx]):(startK-1),];
    
    L1 = matrix(0, nrow= K1*maxLinks,ncol = 2);

    arrowFrom = rep(0,K1*maxLinks);
    arrowTo = rep(0,K1*maxLinks);
    corrTemp = arrowTo;
    patternIndNew = arrowTo; 
    C = CList[[dictIdx]];

    linkInfo = linkDict2Hungary(D1,D2,C,maxLinks = maxLinks,threshold = threshold,template = template);

    p = 1;
    # calculate the coordinates of the second dictionary
    for (j in 1:numRows[dictIdx]){
    	i = patternInd[j];
	indTemp = which(linkInfo[,1] == i);

    	numLinks = length(indTemp); 

    	ixTemp = p:(p+numLinks-1);     
	L1[ixTemp,] = linkDictCoordSingle(L0[j,],numLinks,gap = gaps[dictIdx]);
	
	patternIndNew[ixTemp] = linkInfo[indTemp,2];
	p = p+numLinks;
    }
    
    patternIndNew = patternIndNew[1:(p-1)];
    numRows[dictIdx+1] = p-1;
    L1 = L1[1:(p-1),];
    L = rbind(L,L1);
    patternIndList = c(patternIndList, list(patternIndNew));
    p = 1;
    # copy the images from the second dictionary
    for (I in 1:numRows[dictIdx]){
    	i = patternInd[I];
	indTemp = which(linkInfo[,1] == i);
	indNonZeros = linkInfo[indTemp,2];	
	numLinks = length(indNonZeros);
	
	Y[ind,startK-1+(1:numLinks)] = D2[,indNonZeros];

	

	if (dictIdx == 1){
	   sumTemp = I;
	   }else{
	   sumTemp = I + sum(numRows[1:(dictIdx-1)]);
	   }
    	arrowFrom[p-1+(1:numLinks)] = sumTemp;
    	arrowTo[p-1+(1:numLinks)] = startK-1+(1:numLinks);
    	corrTemp[p-1+(1:numLinks)] = C[i,indNonZeros];
    
	for (j in 1:numLinks){
    	    imgNames[startK+j-1] = paste('p',indNonZeros[j],sep='');
	}
    	startK = startK + numLinks;
    	p = p+numLinks;
    }
    arrowFromList = c(arrowFromList,list(arrowFrom));
    arrowToList = c(arrowToList,list(arrowTo));
    corrTempList = c(corrTempList, list(corrTemp));
}

L[,1] = L[,1]*3;
L[,2] = L[,2]*6;
Y = Y[,1:(startK-1)];
imgNames = imgNames[1:(startK-1)];

out = plotImgOn2dPlane(L,Y,width=32,height=16,scale1=30,scale2=.8,margin = 30);

png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
text(L[,1],L[,2],imgNames,cex = 1,col='red',pos = 1);

for (dictIdx in 1:(numDict-1)){
    arrowFrom = arrowFromList[[dictIdx]];
    arrowTo = arrowToList[[dictIdx]];
    corrTemp = corrTempList[[dictIdx]];
    
    for (i in 1:length(arrowFrom)){
    	x0 = L[arrowFrom[i],1]+ width/(2*dim(out$img)[2]);
    	y0 = L[arrowFrom[i],2];
    	x1 = L[arrowTo[i],1] - width/(2*dim(out$img)[2]);
    	y1 = L[arrowTo[i],2];
    	if (arrowTo[i]*arrowFrom[i]==0){break;}
    	   arrows(x0,y0,x1,y1);
    	   text((x0+x1)/2,(y0+y1)/2,signif(corrTemp[i],2),col='red');
    }
}

if (!is.null(coefVal)){
    
    patternInd = patternIndList[[1]];
    numLessPat = numPatterns;
    for (i in 1:length(numPatterns)){
    	numLessPat[i] = length(patternIndList[[i]]);
    } 

    for (i in 2:length(numPatterns)){
    	patternInd = c(patternInd,patternIndList[[i]]+sum(numPatterns[1:(i-1)]));
    }

    coefTemp = rep(0,sum(numLessPat));
    for (i in 1:length(coefVal)){
    	 temp = which(patternInd==i);
	 if (length(temp)>0){
	    coefTemp[temp] = coefVal[i];
	    }
    }     
    posInd = which(coefTemp>0);
    negInd = which(coefTemp<0);
    col = rep('#000000', length(coefTemp));
    col[posInd] = '#0000FF50';
    col[negInd] = '#FF000050';
    cex = 20*abs(coefTemp);
    points(L[,1],L[,2],col = col,cex = cex,pch = 16);
}



dev.off();



}





linkDict<-function(D1,D2,C,maxLinks = 3,threshold = 0.5,width = 32, height =16,template = template, noValue = FALSE,noNum = FALSE,savePath = NULL){

ind = which(template[,,1]==1);

K1 = dim(D1)[2];
K2 = dim(D2)[2];

Y = matrix(1/3,nrow = width*height,ncol=K1+K1*maxLinks);
imgNames = rep(' ', K1+K1*maxLinks);
L = linkDictCoord(K1,maxLinks=maxLinks);
L0 = L;
L0[,1] = L[,1]*3;
L0[,2] = L[,2]*2;
xCoord = L0[,1];
yCoord = L0[,2];

# fill in the first column with D1
k = 1;
for (i in 1:K1){
    Y[ind,i] = D1[,i];
    imgNames[k] = paste('p',i,sep='');
    k = k + 1;
}
startK = k;
arrowFrom = rep(0,K1*maxLinks);
arrowTo = rep(0,K1*maxLinks);
corrTemp = arrowTo;
p = 1;
for (i in 1:K1){
    cTemp = C[i,];
    temp = sort.int(cTemp,decreasing =TRUE, index.return = TRUE);
    indTemp = temp$ix;
    valTemp = temp$x;
    maxIndex = min(which(valTemp<=threshold))-1; 
    for (j in 1:min(maxIndex,maxLinks)){
    	#print(startK);
	#print(indTemp[j]);
    	Y[ind,startK+j-1] = D2[,indTemp[j]];
	imgNames[startK+j-1] = paste('p',indTemp[j],sep='');
	arrowFrom[p] = i;
	arrowTo[p] = startK+j-1;
	corrTemp[p] = valTemp[j];
	p = p+1;
    }
    startK = startK + maxLinks;
}

     	
out = plotImgOn2dPlane(L0,Y,width=32,height=16,scale1=30,scale2=1);
png(file = savePath,height = dim(out$img)[1]*2,width = dim(out$img)[2]*2);
imagesc(out$img,invert = FALSE);
L = out$L;
text(L[,1],L[,2],imgNames,cex = 1,col='red',pos = 1);

for (i in 1:length(arrowFrom)){
    x0 = L[arrowFrom[i],1]+ width/(2*dim(out$img)[2]);
    y0 = L[arrowFrom[i],2];
    x1 = L[arrowTo[i],1] - width/(2*dim(out$img)[2]);
    y1 = L[arrowTo[i],2];
    if (arrowTo[i]*arrowFrom[i]==0){break;}
    arrows(x0,y0,x1,y1);
    text((x0+x1)/2,(y0+y1)/2,signif(corrTemp[i],2),col='red');
    
}

dev.off();

}


