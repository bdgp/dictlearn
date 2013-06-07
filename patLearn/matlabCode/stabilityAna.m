% The following matlab code attempts to perform stability analysis
% of the dictionary learning algorithm with respective to various
% perturbation: e.g. different initial values.

addpath('~/matlab/packages/spams-matlab/');
addpath('~/matlab/packages/spams-matlab/build/');
addpath('~/fruitfly/image/osDict/utilities/');
addpath('../CV');
addpath('../');

numReplicates = 100;

sigma = 0;
numPatterns = 1:2:100;

resolution = '32by16';
if strmatch(resolution,'64by32')
    initial_path = ['./64by32/randomStartSigma=',num2str(sigma),'/'];
    width = 64;
    height = 32;
else
    initial_path = ['./64by32/randomStartSigma=',num2str(sigma),'/'];
    width = 32;
    height = 16;
end


doObj = 1 % plot the objective function for all sizes and repetitions?
if doObj
    objFcnVal = zeros(numReplicates,length(numPatterns));
    for k = 1:length(numPatterns)
        K = numPatterns(k);
        current_path0 = [initial_path, 'K=', num2str(K),'/']; 
        load([current_path0,'bestDict.mat']);
        objFcnVal(:,k) = R';
    end
    figure; plot(1:numReplicates,objFcnVal);
    xlabel('replicates');
    ylabel('objective function value');
    print(gcf,'-dpng',[initial_path,'objFcnVal1.png']);

    figure; plot(numPatterns,objFcnVal');
    xlabel('K');
    ylabel('objective function value');
    print(gcf,'-dpng',[initial_path,'objFcnVal2.png']);
    
end

    
loadResults = 1;
if loadResults
    estStability = zeros(1,length(numPatterns));
    estDictError = estStability;

for k = 1:length(numPatterns)
    K = numPatterns(k);

    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    load([current_path0,'bestDict.mat']);
    D0 = Dbest;
    
    Dhat = zeros(size(D0,1),K,numReplicates);

    for L = 1:numReplicates
        load([current_path0,'rep',num2str(L), ...
                      'Dict.mat']);
        Dhat(:,:,L) = D;
    end
    
    distMat = zeros(numReplicates,numReplicates);
    errorMat = zeros(numReplicates,numReplicates);
    for q = 1:numReplicates
        for p = q:numReplicates
            [distMat(q,p),distMat(p,q),errorMat(q,p),errorMat(p,q)] = distDictRegression(Dhat(:,:,q),Dhat(:,:,p));
        end
    end
    save([current_path0,'distMatrixDict.mat'],'distMat','errorMat');    
    estStability(k) = sum(distMat(:))/(numReplicates*(numReplicates-1));
    estDictError(k) = sum(errorMat(:))/(numReplicates*(numReplicates-1));

end
save([initial_path,'estStabDict.mat'],'numPatterns','estStability','estDictError');
figure; plot(numPatterns,estStability);
print(gcf,'-dpng',[initial_path,'Dstability.png']);
close all;

figure; plot(numPatterns,estDictError);
print(gcf,'-dpng',[initial_path,'Derror.png']);
close all;

figure; scatter(estStability,estDictError);
temp = num2strBatch(numPatterns);
text(estStability,estDictError, temp);
print(gcf,'-dpng',[initial_path,'stabVSError.png']);

end


doShiftCV = 1;
if doShiftCV
    [m,n] = size(X);
    param.mode = 2; 
    param.lambda=0; 
    param.numThreads=-1; 
    param.batchsize=min(1024,n);
    param.pos = 1; 
    param.posAlpha = 1;
   
    Xpert = X;
    for ii = 1:size(X,2)
        method = randsample(1:4,1);
        pp.by = 2;
        Ytemp = perturbImg(Y(:,ii),width,height,method,pp);
        Xpert(:,ii) = Ytemp(ind,:);
    end

    Lambda = [0,0.1,0.2,0.5,1,2,5];
    cvError = zeros(length(Lambda),length(numPatterns));
    for l = 1:length(Lambda)
        param.lambda = Lambda(l);    
        for k = 1:length(numPatterns)
            K = numPatterns(k);    
            current_path0 = [initial_path, 'K=', num2str(K),'/']; 
            load([current_path0,'bestDict.mat']);
            Dstd = Dbest;
            for q = 1:K
                Dstd(:,q) = Dbest(:,q)/max(Dbest(:,q));
            end
            alpha = mexLasso(X,Dstd,param);
            
            E = Xpert - Dstd*alpha;
            cvError(l,k) = mean(E(:).^2);
        end
    end
end

% find interval estimates of the distMat and errorMat measures
loadError = 0;
if loadError
    estStability = zeros(1,length(numPatterns));
    estDictError = estStability;
    estDictErrorU = estStability;
    estDictErrorL = estStability;
    estStabilityU = estStability;
    estStabilityL = estStability;
    

    for k = 1:length(numPatterns)
        K = numPatterns(k);
        current_path0 = [initial_path, 'K=', num2str(K),'/']; 
        load([current_path0,'distMatrixDict.mat']);
        
        estStability(k) = median(distMat(:));
        estStabilityU(k) = quantile(distMat(:),.75);
        estStabilityL(k) = quantile(distMat(:),.25);
   
        estDictError(k) = median(errorMat(:));
        estDictErrorU(k) = quantile(errorMat(:),.75);
        estDictErrorL(k) = quantile(errorMat(:),.25);
    end

    figure; plot(numPatterns,estStability);
    hold on;
    plot(numPatterns,estStabilityU,'r');
    plot(numPatterns,estStabilityL,'r');
    hold off
    figure; plot(numPatterns,estDictError);
    hold on;
    plot(numPatterns,estDictErrorU,'r');
    plot(numPatterns,estDictErrorL,'r');
    hold off
end

    