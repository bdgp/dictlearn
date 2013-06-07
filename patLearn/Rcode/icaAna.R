res = fastICA(t(X),n.comp=70,maxit = 500);
D = t(res$A);
alpha = t(res$S); # notice that alpha%*%t(alpha) = I
imageBatchDisplaySaveMemory(D,template=template,colorScale='blueRed');
imageBatchDisplaySaveMemory(D%*%alpha,template=template,colorScale='blueRed');