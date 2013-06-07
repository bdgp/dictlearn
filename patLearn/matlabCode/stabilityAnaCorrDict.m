addpath('~/matlab/packages/spams-matlab/');
addpath('~/matlab/packages/spams-matlab/build/');
addpath('~/fruitfly/image/image2TI');
addpath('~/fruitfly/image/Visualization');
addpath('~/fruitfly/image/EDA/explore');
addpath('~/fruitfly/annotation/TF/som/somtoolbox/');
addpath('~/topoDict/');
addpath('~/fruitfly/image/osDict/utilities/');
addpath('../CV');
addpath('../');

%rng(215);
numReplicates = 10;

sigma = 0;
numPatterns = 1:2:100;

initial_path = ['./32by16/randomStartSigma=',num2str(sigma),'/'];


param.mode = 2; 
param.lambda=0; % penalizing constant on alpha
param.numThreads=-1; % number of threads
param.pos = 1; % positive alpha
param.posAlpha = 1;


width = 32;
height = 16;
loadCV = 1;
if loadCV
    load('./32by16/data.mat');
    estStability = zeros(1,length(numPatterns));

for k = 1:length(numPatterns)
    K = numPatterns(k);
    K = 93
    
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    load([current_path0,'bestDict.mat']);
    D0 = Dbest;
    
    Dhat = zeros(size(X,1),K,numReplicates);
    Alpha = zeros(K,size(X,2),numReplicates);
    
    for L = 1:numReplicates
        load([current_path0,'rep',num2str(L), ...
                      'Dict.mat']);
        % should standardize ...
        Dhat(:,:,L) = D;
        Alpha(:,:,L) = mexLasso(X,D,param);
    end
    
    
    distMat = zeros(numReplicates,numReplicates);
    errorMat = zeros(numReplicates,numReplicates);
    for q = 1:(numReplicates-1)
        for p = (q+1):numReplicates
            %Ctemp = corr(Dhat(:,:,q),Dhat(:,:,p));
            %Ctemp = corr(Alpha(:,:,q)',Alpha(:,:,p)');
            [distMat(q,p),distMat(p,q),errorMat(q,p),errorMat(p,q)] = distDictRegression(Alpha(:,:,q)',Alpha(:,:,p)');

            %Ctemp(Ctemp<0.5) = 0;
            % can have multiple strong correlations. local picture only
            %errorTemp = amariError(Ctemp);
            %distMat(q,p) = errorTemp;
            %distMat(p,q) = errorTemp;
        end
    end
    stop
    save([current_path0,'distMatrixAlphaCorr.mat'],'distMat');
    estStability(k) = sum(distMat(:))/(numReplicates*(numReplicates-1));    

end
save([initial_path,'estStabAlphaCorrDict.mat'],'numPatterns','estStability');
figure; plot(numPatterns,estStability);
print(gcf,'-dpng',[initial_path,'DstabilityCorrAlpha.png']);
close all;


end

