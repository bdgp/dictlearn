addpath('~/matlab/packages/spams-matlab/');
addpath('~/matlab/packages/spams-matlab/build/');
addpath('~/fruitfly/image/image2TI');
addpath('~/fruitfly/image/Visualization');
addpath('~/fruitfly/image/EDA/explore');
addpath('~/fruitfly/annotation/TF/som/somtoolbox/');
addpath('~/topoDict/');
addpath('~/fruitfly/image/osDict/utilities/');
addpath('../');
addpath('../CV/');

if ~exist('Load')
    resolution = '32by16';
    if strmatch(resolution,'64by32')
        load('../dataTemp.mat');
        width = 64;
        height = 32;
        template = generateTemplate();
        template = imresize(template,.5,'nearest');
    else
        load('./32by16/data.mat');
        width = 32;
        height = 16;
    end
    load('./geneNames.mat');
    load('~/fruitfly/annotation/TF/annotDataFinal.mat');
    
    tfNames = gsym;
    tfInd = [];
    for i = 1:length(tfNames)
        indTemp = strmatch(tfNames(i),geneName,'exact');
        tfInd = [tfInd;indTemp];
    end
    ind = find(template(:,:,1)==1);
    %X = X(:,tfInd);

    Y = zeros(height*width,size(X,2));
    Y(ind,:) = X;
    
    Load = 1;
end

sigma = 0;
initial_path = ['./32by16/randomStartSigma=',num2str(sigma),'/'];

[m,n] = size(X);
param.mode = 2; 
param.lambda=0; % penalizing constant on alpha
param.numThreads=-1; % number of threads
param.pos = 1; % positive alpha
param.posAlpha = 1;
numPixels = length(ind);
rep = 10;
numPatterns = 1:2:100;
cvError = zeros(length(numPatterns),rep);
fold = 2;
mskTemp = squareMask(width,height,8);
msk = mskTemp(:);
msk = msk(ind);
testInd = find(msk == 1);
trainInd = setdiff(1:length(msk),testInd);
indSet{1} = testInd;
indSet{2} = trainInd;

for k = 1:length(numPatterns)
    K = numPatterns(k)
    testError = zeros(rep,1);
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    %load([current_path0,'bestDict.mat']);
    %D = Dbest;
    
    for i = 1:rep
        load([current_path0,'rep',num2str(i),'Dict.mat']);   
        testErrorTemp = zeros(fold,n);
        for f = 1:fold
                testInd = indSet{f};
                trainInd = setdiff(1:numPixels,testInd);
                alphaTrain = mexLasso(X(trainInd,:),D(trainInd,:),param);
                E = X(testInd,:) - D(testInd,:)*alphaTrain;
                testErrorTemp(f,:) = mean(E.^2);
            
        end
        testError(i) = mean(testErrorTemp(:));
        
    end
    cvError(k,:) = testError';
end

figure;
boxplot(cvError',numPatterns);
xlabel('K');
ylabel('test error');
title('pixel deletion: 2 by 2 square block');
stop

for k = 1:length(numPatterns)
    K = numPatterns(k)
    testError = zeros(rep,n);
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    load([current_path0,'bestDict.mat']);
    D = Dbest;
    rng('default')
    for i = 1:rep
        testErrorTemp = zeros(fold,n);
        for imgIdx = 1:n
        
            indSet = testTrainIndProduce(numPixels,fold,'random');
            for f = 1:fold
                testInd = indSet{f};
                trainInd = setdiff(1:numPixels,testInd);
                alphaTrain = mexLasso(X(trainInd,imgIdx),D(trainInd,:),param);
                E = X(testInd,imgIdx) - D(testInd,:)*alphaTrain;
                testErrorTemp(f,imgIdx) = mean(E.^2);
            end
        end
        testError(i,:) = sum(testErrorTemp,1);
        
    end
    cvError(k,:) = mean(testError,1);
end


stop



for k = 1:length(numPatterns)
    K = numPatterns(k);
    testError = zeros(rep,n);
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    load([current_path0,'bestDict.mat']);
    D = Dbest;
    rng('default')
    for i = 1:rep
        trainInd = randsample(numPixels,floor(numPixels/2));
        testInd = setdiff(1:numPixels,trainInd);
        alphaTrain = mexLasso(X(trainInd,:),D(trainInd,:),param);
        E = X(testInd,:) - D(testInd,:)*alphaTrain;
        E1 = mean(E.^2);
        alphaTest = mexLasso(X(testInd,:),D(testInd,:),param);
        E = X(trainInd,:) - D(trainInd,:)*alphaTest;
        E2 = mean(E.^2);
        testError(i,:) = E1 + E2;
      
        
    end
    cvError(k,:) = mean(testError);
end


stop

numPixels = length(ind);
rep = 10;
Lambda = [0.0005,0.001,0.002,0.005,0.01,0.02,0.05];
cvError = zeros(length(Lambda),n);

for j = 1:length(Lambda)
    lambda = Lambda(j);
    param.lambda = lambda;
    testError = zeros(rep,n);
    
    for i = 1:rep
        trainInd = randsample(numPixels,floor(numPixels/2));
        testInd = setdiff(1:numPixels,trainInd);
        alphaTrain = mexLasso(X(trainInd,:),D(trainInd,:),param);
        E = X(testInd,:) - D(testInd,:)*alphaTrain;
        testError(i,:) = mean(E.^2);
    end
    cvError(j,:) = mean(testError);
end