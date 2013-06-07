addpath('~/matlab/packages/spams-matlab/');
addpath('~/matlab/packages/spams-matlab/build/');
addpath('~/fruitfly/image/image2TI');
addpath('~/fruitfly/image/Visualization');
addpath('~/fruitfly/image/EDA/explore');
addpath('~/fruitfly/annotation/TF/som/somtoolbox/');
addpath('~/topoDict/');
addpath('~/fruitfly/image/osDict/utilities/');
addpath('../CV');


if ~exist('Load')
load('./32by16/data.mat');
ind = find(template(:,:,1)==1);
Load = 1;
end

width = 32; height = 16;
iter = 300;
numPatterns = 1:2:100;
numReplicates = 2;
n = size(X,2);
initial_path = ['./32by16/CV/',num2str(numReplicates),'fold/'];


for k = 1:length(numPatterns)
    K = numPatterns(k);
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    mkdir(current_path0);
    
    pathTemp = ['./32by16/randomStart/K=', num2str(K),'/']; 
    load([pathTemp,'bestDict.mat']);
    D0 = Dbest;
    lambda = 0;
    gamma1 = 0;

    [m,n] = size(X);
    param.mode = 2; 
    param.K=K;  % learns a dictionary with K elements
    param.lambda=lambda; % penalizing constant on alpha
    param.numThreads=-1; % number of threads
    param.batchsize=min(1024,n);
    param.posD = true;   % positive dictionary
    param.iter = 300;  % number of iteration; default = 1000
    param.modeD = 0;
    param.verbose = 0;
    param.pos = 1; % positive alpha
    param.posAlpha = 1;
    param.gamma1 = gamma1; % penalizing constant on dict    
    
    param.ols = 1; % this might be wrong!
    
    param.D = D0; 
    %%%%%%%%%% EXPERIMENT %%%%%%%%%%%
  
    % should make deterministic? 
    rng('default');
    indSet = testTrainIndProduce(size(X,2),numReplicates,'kmeans');
    for L = 1:numReplicates
        testInd = indSet{L};
        trainInd = setdiff((1:n)',testInd);
        DTrain = mexTrainDL(X(:,trainInd),param);
        alpha_test = mexLasso(X(:,testInd),DTrain,param);
        testError=mean(0.5*sum((X(:,testInd)-DTrain*alpha_test).^2));

                  
                  %mexLasso(X(:,trainInd),DTrain,param);
                  %idx = randsample(K,K);
                  %DTrain = DTrain(:,idx);
                  %C = dissimilarityDict(DTrain,DTest,'euclidean');
                  %[ass,cost_euclidean(L)] = munkres(C);
                 
                  D = DTrain;
                  Dtemp = zeros(width*height,K) + max(D(:))/4;      
                  Dstd = dictStd(D,2);
                  Dtemp(ind,:) = Dstd;
                  
                  imgMat = imageBatchProduce(Dtemp,width,height,floor(sqrt(K))+1,floor(sqrt(K))+1);
                  M = max(imgMat(:));
                  if M > 0
                      scale = 1;
                      imgMat(imgMat>=M/scale) = M/scale;
                      imgMat = imgMat/(M/scale);
                  end
                  imwrite(1-imgMat,[current_path0,'cvRep',num2str(L),'.jpg'],'jpg');
     
                  %{
                  D1 = zeros(width*height,K) + max(DTrain(:))/4;      
                  Dstd1 = dictStd(DTrain,2);
                  D1(ind,:) = Dstd1;

                  D2 = zeros(width*height,K) + max(DTest(:))/4;      
                  Dstd2 = dictStd(DTest,2);
                  D2(ind,:) = Dstd2;
                  imageBatchDisplay(D1, width, height, 5,5,'');
                  imageBatchDisplay(D2(:,ass), width, height, 5,5,'');
                  
                  
                  
                  input('press enter to continue');
                  close all;
                  %}
                  param.D = [];
                  save([current_path0,'cvRep',num2str(L),'.mat'], 'DTrain','param','trainInd','testInd','testError');     
                  %set(gcf,'PaperPositionMode','auto'); % preserve the size
                  %print(gcf,'-dpng',[current_path3,'spc',num2str(L),'.png']);
      
                  %close(gcf);
                  param.D = D0;
                  
    end
              
    input('press enter to continue');
              
end
