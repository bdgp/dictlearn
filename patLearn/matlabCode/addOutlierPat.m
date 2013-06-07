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

Y = crossPat(width,height,6,1);
Y = Y(:);
Y = Y(ind);
numOut = 2;
Xextra = repmat(Y,1,numOut);
[m,n] = size(X);
lambda = 0;
gamma1 = 0;
numPatterns = 1:4:100;
replicates = 10;
corrTemp = zeros(length(numPatterns),replicates);
%{
for k = 1:length(numPatterns)
    K = numPatterns(k);

    path = ['./',resolution,'/addOutlierPat/numOutlier=',num2str(numOut),'/K=',num2str(K),'/'];
    mkdir(path);
   
    
   
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
    D0 = dictLearnInit(X,K,'random',0);
    param.D = D0; 
         
    Dtemplate = mexTrainDL([X,Xextra],param);
    %D2 = dictLearn(X,param);
    %stop
    maxCorr = zeros(1,replicates);
    for i = 1:replicates
        D0 = dictLearnInit(X,K,'random',0);
        param.D = D0; 
        %%%%%%%%%% EXPERIMENT %%%%%%%%%%%
        tic
        D = mexTrainDL([X,Xextra],param); 
        
        permInd = dictMatch(Dtemplate,D,'euclidean');
      
        toc
        D = D(:,permInd);
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
        imwrite(1-imgMat,[path,'rep',num2str(i),'.jpg']);
      
        save([path,'rep',num2str(i),'Dict.mat'],'D');

        input('press enter to continue');

        %imageBatchDisplaySaveMemory(D,32,16,floor(sqrt(K))+1,floor(sqrt(K))+1,'',.25);
        maxCorr(i) = max(corr(Y,D));
    end 
    save([path,'maxCorr.mat'],'maxCorr');
    corrTemp(k,:) = maxCorr;
end
%}
corrTemp = zeros(length(numPatterns),replicates);
for k = 1:length(numPatterns)
    K = numPatterns(k);
    path = ['./',resolution,'/addOutlierPat/numOutlier=',num2str(numOut),'/K=',num2str(K),'/'];
    load([path,'maxCorr']);
    corrTemp(k,:) = maxCorr;
end
path2 = ['./',resolution,'/addOutlierPat/numOutlier=',num2str(numOut),'/maxCorr.mat'];
save(path2,'corrTemp');
prob = mean(corrTemp>.9,2);
figure;scatter(numPatterns,prob);
xlabel('K');
ylabel(['probability of occurrence of the artificial pattern in ',num2str(replicates),' runs'])
print(gcf,'-dpng',[path2,'probOccur.png']);


