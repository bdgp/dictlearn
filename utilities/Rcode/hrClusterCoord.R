hrClusterCoord<-function(hr){
height = hr$height;
merge = hr$merge;
mergeOrder = rep(0,length = (dim(merge)[1])+1);
mergeHeight = mergeOrder;
k = 1;
for (i in 1:(dim(merge)[1])){
    #if (k==51){print(mergeOrder);}
    if (merge[i,1]<0){
       mergeOrder[k] = -merge[i,1];
       mergeHeight[k] = height[i];
       k = k + 1;
    }
    if (merge[i,2]<0){
       mergeOrder[k] = -merge[i,2];
       mergeHeight[k] = height[i];
       k = k + 1;
    }
 }
temp1 = mergeOrder;
yCoord = xCoord = mergeOrder;

for (i in 1:length(mergeOrder)){
    temp1[which(hr$order == mergeOrder[i])] = mergeHeight[i];
    xCoord[i] = which(hr$order == i);
    
}

mergeHeight[hr$order] = temp1;
yCoord = mergeHeight;
L = cbind(xCoord,yCoord);
return(L);
}