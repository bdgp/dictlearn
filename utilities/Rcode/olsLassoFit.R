olsLassoFit <- function(X,D,alpha,intercept = TRUE){
n = dim(X)[2];
if (intercept){
   alphaOLS = rbind(alpha,rep(1,dim(alpha)[2]));
   Dnew = cbind(D,rep(1,dim(X)[1]));
}else{
   alphaOLS = alpha;
   Dnew = D;
}		

rSq = rep(0,n);
for (i in 1:n){
    res = olsLassoFit0(X[,i],Dnew,alphaOLS[,i]);
    alphaOLS[,i] = res$alphaOLS;
    rSq[i] = res$r.squared;
}

return(list(alphaOLS = alphaOLS,rSq = rSq, D = Dnew));
}

olsLassoFit0<-function(x,D,a){
nonZeroInd = which(abs(a)>1e-6);
Dnew = D[,nonZeroInd];
lm.res = lm(x~Dnew-1);
s = summary(lm.res);
r.squared = s$r.squared;
coef = lm.res$coefficients;
alphaOLS = a;
alphaOLS[nonZeroInd] =coef;
return(list(alphaOLS=alphaOLS,r.squared = r.squared));
}