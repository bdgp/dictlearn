drawDendrogram<-function(hr,mergeHeight,L){
merge = hr$merge;
plotOrder = hr$order;
xCoord = sort(L[,1]);

clusterCoords = matrix(0,nrow = length(mergeHeight),ncol = 2);
for (i in 1:length(mergeHeight)){
     #print(i)
     if (merge[i,1]<0 & merge[i,2]<0){
     	whichPt1 = -merge[i,1];
     	whichPt2 = -merge[i,2];
     	temp1 = xCoord[which(plotOrder == whichPt1)];
     	temp2 = xCoord[which(plotOrder == whichPt2)];
	
     	x = (temp1 + temp2)/2;
     	y = mergeHeight[i];
     	clusterCoords[i,] = c(x,y);
	segments(temp1,y,temp2,y,col='blue',lwd =2);		  
	
     }else if (merge[i,1]<0 & merge[i,2]>0){
     	whichPt1 = -merge[i,1];
     	whichCluster = merge[i,2];

     	temp1 = xCoord[which(plotOrder == whichPt1)];
     	temp2 = clusterCoords[whichCluster,1];
	
     	x = (temp1 + temp2)/2;
     	y = mergeHeight[i];
     	clusterCoords[i,] = c(x,y);
	segments(temp1,y,temp2,y,col='blue');
	
	clusterHeight = clusterCoords[whichCluster,2];
	if (clusterHeight<y){
	   segments(temp2,clusterHeight,temp2,y,col='blue');
	
	}else{
	   segments(temp1,clusterHeight,temp1,y,col='blue');
	}		  
	
     }else if (merge[i,1]>0 & merge[i,2]<0){
     	
     	whichCluster = merge[i,1];
     	whichPt2 = -merge[i,2];
     	temp2 = xCoord[which(plotOrder == whichPt2)];
     	temp1 = clusterCoords[whichCluster,1];
     	x = (temp1 + temp2)/2;
     	y = mergeHeight[i];
     	clusterCoords[i,] = c(x,y);

     	segments(temp1,y,temp2,y,col='blue');
	
	clusterHeight = clusterCoords[whichCluster,2];
	if (clusterHeight<y){
	   segments(temp1,clusterHeight,temp1,y,col='blue');
	
	}else{
	   segments(temp2,clusterHeight,temp2,y,col='blue');
	}		  
		  
	}else if (merge[i,1]>0 & merge[i,2]>0){
	
     	whichCluster1 = merge[i,1];
     	whichCluster2 = merge[i,2];
       	temp1 = clusterCoords[whichCluster1,1];
    	temp2 = clusterCoords[whichCluster2,1];
     	x = (temp1 + temp2)/2;
     	y = mergeHeight[i];
     	clusterCoords[i,] = c(x,y);
	segments(temp1,y,temp2,y,col='blue');		  
	
	clusterHeight1 = clusterCoords[whichCluster1,2];
	clusterHeight2 = clusterCoords[whichCluster2,2];
	
	segments(temp2,clusterHeight2,temp2,y,col='blue');
	segments(temp1,clusterHeight1,temp1,y,col='blue');
	
	
     }
   }
}