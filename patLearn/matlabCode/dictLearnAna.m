% This is the matlab code to start with...

addpath('../spams-matlab/');
addpath('../spams-matlab/build/');
addpath('../../utilities/');

if ~exist('Load')
    % choose the size of the dictionary
    % either 64 by 32 or 32 by 16
    
    resolution = '32by16';
    if strmatch(resolution,'64by32')
        load('../../data/data64by32.mat');
        width = 64;
        height = 32;
        ratio = .5;
    else
        load('../../data/data32by16.mat');
        width = 32;
        height = 16;
        ratio = .25;
    end
    
    readTF = 0; % read the TF annotation data?
    if readTF
        load('../../data/annotDataFinal.mat');
        tfNames = gsym;
        tfInd = [];
        for i = 1:length(tfNames)
            indTemp = strmatch(tfNames(i),geneNames,'exact');
            tfInd = [tfInd;indTemp];
        end
    end
    ind = find(template(:,:,1)==1);
    Y = zeros(height*width,size(X,2));
    Y(ind,:) = X;
    Load = 1;
end

[m,n] = size(X);

K = 49;

% regularization parameters
lambda = 0.1; % sparsity on the coefficients alpha
gamma1 = 0.1; % sparsity on the dictionary patterns

% please refer to Julien's documentation for the package spams...      
param.mode = 2;
param.modeD = 1; 

param.K=K;  
param.lambda=lambda; 
param.numThreads=-1; 
param.batchsize=min(1024,n);
param.posD = true;   % positive dictionary
param.iter=300;     
param.verbose = 1; % print out the iteration
param.pos = 1; % positive alpha
param.posAlpha = 1;
param.gamma1 = gamma1; % penalizing constant on dict
      
% initial value: randomly select K images from the data set
D0 = dictLearnInit(X,K,'random',0);
param.D = D0; 
     
%tic              
D = mexTrainDL(X,param); 
alpha = mexLasso(X,D,param);
%toc              
imageBatchDisplaySaveMemory(D,width,height,floor(sqrt(K))+1,floor(sqrt(K))+1,'',ratio);
      

