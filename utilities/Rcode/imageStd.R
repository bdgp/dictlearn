imageStd <- function(X,method = 'max',lower = 0.5){
# Standardization of expression images
# X: data
# method: 
# 1) max: x/max(x);
# 2) mean: (x-mean(x))_+;
# 3) quantile: (x - quantile(x,lower))_+;
# 4) L1: x/||x||_1;
# 5) L2: x/||x||_2;

Xout = X;
if (method == 'max'){
   for (i in 1:dim(X)[2]){
     Xout[,i] = X[,i]/max(X[,i]);
   }
}else if (method == 'mean'){
   for (i in 1:dim(X)[2]){
      temp = (X[,i] - mean(X[,i])); temp[temp<0]=0;
      Xout[,i] = temp;
    }
}else if (method == 'quantile'){
   for (i in 1:dim(X)[2]){
      temp = (X[,i] - quantile(X[,i],lower)); temp[temp<0]=0;    
      Xout[,i] = temp;
   }
}else if (method == 'L1'){
   for (i in 1:dim(X)[2]){
       Xout[,i] = X[,i]/sum(abs(X[,i]));
   }
}else if (method == 'L2'){
   for (i in 1:dim(X)[2]){
       Xout[,i] = X[,i]/sqrt(sum(X[,i]^2));
}
}else{return(NULL);}
return(Xout);     
   

}    
