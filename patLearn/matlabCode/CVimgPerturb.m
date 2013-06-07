addpath('~/matlab/packages/spams-matlab/');
addpath('~/matlab/packages/spams-matlab/build/');
addpath('~/fruitfly/image/image2TI');
addpath('~/fruitfly/image/Visualization');
addpath('~/fruitfly/image/EDA/explore');
addpath('~/fruitfly/annotation/TF/som/somtoolbox/');
addpath('~/topoDict/');
addpath('~/fruitfly/image/osDict/utilities/');
addpath('../');
%iter = [200,400,600,800,1000]; % K = 25, 36, 49, 64, 81

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

numPatterns = 1:2:100;
[m,n] = size(X);
param.mode = 2; 
param.lambda=0; % penalizing constant on alpha
param.numThreads=-1; % number of threads
param.pos = 1; % positive alpha
param.posAlpha = 1;


  
numReplicates = 10;
XpertEm = zeros(size(X,1),size(X,2),numReplicates);
Lambda = [0,0.1,0.2,0.5,1,2,5];
cvError = zeros(length(numReplicates),length(numPatterns));
Prob = 0.1:0.1:1
for pp = 1:length(Prob)
    prob = Prob(pp);
for L = 1:(numReplicates)
    for ii = 1:size(X,2)
        change = rand(1);
        if change < prob
            method = randsample(1:4,1);
            pp.by = 1;
            Ytemp = perturbImg(Y(:,ii),width,height,method,pp);
            XpertEm(:,ii,L) = Ytemp(ind,:);
        else
            XpertEm(:,ii,L) = X(:,ii);
        end
    end
end

% ESCV
for k = 1:length(numPatterns)
        K = numPatterns(k)

        current_path0 = [initial_path, 'K=', num2str(K),'/']; 
        load([current_path0,'bestDict.mat']);
        D = Dbest;
        Dstd = D;
        for q = 1:K
            Dstd(:,q) = D(:,q)/max(D(:,q));
        end
        
        Y1 = zeros(size(X,1),size(X,2),numReplicates);
        for L = 1:numReplicates            
            Xtemp = XpertEm(:,:,L);
            alpha = mexLasso(Xtemp,Dstd,param);
            Y1(:,:,L) = Dstd*alpha;
        end
        Ymean = mean(Y1,3);
        Ysq = mean(Y1.^2,3);
        escvVar(k) = sum(sum(Ysq - Ymean.^2));
        escvError(k) = escvVar(k)/sum(Ysq(:));
end
    path0 = './32by16/CV/imgPerbCV/KSelectionESCV/';
    mkdir(path0);
    save([path0,num2str(prob),'Shift.mat'],'escvError','param', ...
         'prob','escvVar');
    figure;plot(numPatterns,escvError);
    print(gcf,'-dpng',[path0,num2str(prob),'Shift.png']);
    figure;plot(numPatterns(20:end),escvError(:,20:end));
    print(gcf,'-dpng',[path0,num2str(prob),'ShiftZoom.png']);
    close all;

end

%{
%for l = 1:length(Lambda)
%    param.lambda = Lambda(l);
    %param.lambda = .2;
    for k = 1:length(numPatterns)
        K = numPatterns(k)
        %K = 59;
        current_path0 = [initial_path, 'K=', num2str(K),'/']; 
        load([current_path0,'bestDict.mat']);
        D = Dbest;
        Dstd = D;
        for q = 1:K
            Dstd(:,q) = D(:,q)/max(D(:,q));
        end
        
        for L = 1:numReplicates            
            Xtemp = XpertEm(:,:,L);
            alpha = mexLasso(Xtemp,Dstd,param);
            testSampleInd = setdiff(1:numReplicates,L);
            errorTemp = 0;
            for j = 1:(numReplicates-1)
                Xpert = XpertEm(:,:,testSampleInd(j));
                E = Xpert - Dstd*alpha;
                errorTemp = errorTemp + mean(E(:).^2);
            end
            
            cvError(L,k) = errorTemp/numReplicates;
        end
    end
    path0 = './32by16/CV/imgPerbCV/KSelection/';
    
    save([path0,num2str(prob),'Shift.mat'],'cvError','param', ...
         'prob');
    figure;boxplot(cvError,numPatterns);
    print(gcf,'-dpng',[path0,num2str(prob),'Shift.png']);
    figure;boxplot(cvError(:,20:end),numPatterns(20:end));
    print(gcf,'-dpng',[path0,num2str(prob),'ShiftZoom.png']);
    close all;
    %}

%{
% load multiresolution dictionary here...
loadDict = true;
if loadDict
    sigma = 0
    initial_path = ['./32by16/randomStartSigma=',num2str(sigma),'/'];
    initial_path = ['./32by16/randomStart/'];
    numPatterns = [6,13,25,47,1];
    D = zeros(length(ind),sum(numPatterns));
    q = 1;
    for i = 1:length(numPatterns)
        K = numPatterns(i);
        current_path0 = [initial_path, 'K=', num2str(K),'/']; 
        load([current_path0,'bestDict.mat']);
        D(:,q:(q+K-1)) = Dbest;
        q = q+K;
    end
    K = size(D,2);
end

[m,n] = size(X);
param.mode = 2; 
param.K=K;  % learns a dictionary with K elements
param.numThreads=-1; % number of threads
param.batchsize=min(1024,n);
param.verbose = 0;

Lambda = [0.01,0.02,0.05,0.1,0.2,0.5,1,2,5];
estStabAlpha = zeros(n,length(Lambda));
% estimation stability alpha:
% for each perturbed images, estimate the coefficients and look at
% the variability of them.
%{
for i = 1:length(Lambda)
    lambda = Lambda(i);
    param.lambda = lambda;
    alphaMat = zeros(K,n,4);
    
    
    p.by = 1;
    p.direction = 'vertical';
    Xpert1 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert1(ind,:);
    alpha1 = mexLasso(Xtemp,D,param);
    alphaMat(:,:,1) = alpha1;
    
    p.by = -1;
    Xpert2 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert2(ind,:);
    alpha2 = mexLasso(Xtemp,D,param);
    alphaMat(:,:,2) = alpha2;
   
    p.by = 1;
    p.direction = 'horizontal';
    Xpert3 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert3(ind,:);
    alpha3 = mexLasso(Xtemp,D,param);
    alphaMat(:,:,3) = alpha3;

    
    p.by = -1;
    Xpert4 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert4(ind,:);
    alpha4 = mexLasso(Xtemp,D,param);
    alphaMat(:,:,4) = alpha4;

    varAlpha = zeros(n,1);
    sqAlpha = zeros(n,1);
    for j = 1:n
        alphaTemp = zeros(K,4);
        for k = 1:4
            alphaTemp(:,k) = alphaMat(:,j,k);
        end
        covTemp = cov(alphaTemp');
        varAlpha(j) = sum(diag(covTemp));
        sqAlpha(j) = sum(alphaTemp(:).^2);
        
    end
    estStabAlpha(:,i) = varAlpha./sqAlpha;
end
%}
%{
shift = 1;
cvError = zeros(n,length(Lambda));
for i = 1:length(Lambda)
    lambda = Lambda(i);
    param.lambda = lambda;
    alpha = mexLasso(X,D,param);
    Xest = D*alpha;
    
    p.by = shift;
    p.direction = 'vertical';
    Xpert1 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert1(ind,:);
    errorTemp1 = Xest - Xtemp;    
    
    p.by = -shift;
    Xpert2 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert2(ind,:);
    errorTemp2 = Xest - Xtemp;
   
    p.by = shift;
    p.direction = 'horizontal';
    Xpert3 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert3(ind,:);
    errorTemp3 = Xest - Xtemp;

    
    p.by = -shift;
    Xpert4 = perturbImg(Y,width,height,'shift',p);
    Xtemp = Xpert4(ind,:);
    errorTemp4 = Xest - Xtemp;

    cvError(:,i) = (sum(errorTemp1.^2)+sum(errorTemp2.^2)+sum(errorTemp3.^2)+sum(errorTemp4.^2))';
end
figure;boxplot(cvError,Lambda);
xlabel('lambda');
ylabel('CV error');
%}

shift = 2;
cvError = zeros(n,length(Lambda));
Xall = zeros(m,n,5);
Xall(:,:,1) = X;

p.by = shift;
p.direction = 'vertical';
Xpert1 = perturbImg(Y,width,height,'shift',p);
Xall(:,:,2) = Xpert1(ind,:);
    
p.by = -shift;
Xpert2 = perturbImg(Y,width,height,'shift',p);
Xall(:,:,3) = Xpert2(ind,:);
   
p.by = shift;
p.direction = 'horizontal';
Xpert3 = perturbImg(Y,width,height,'shift',p);
Xall(:,:,4) = Xpert3(ind,:);
    
p.by = -shift;
Xpert4 = perturbImg(Y,width,height,'shift',p);
Xall(:,:,5) = Xpert4(ind,:);



for i = 1:length(Lambda)
    lambda = Lambda(i);
    param.lambda = lambda;
    for j = 1:5
        Xtrain = Xall(:,:,j);
        alpha = mexLasso(Xtrain,D,param);
        Xest = D*alpha;
        errorTemp = zeros(m,n,4);
        q = 1
        for k = 1:5
            if k~=j
                errorTemp(:,:,q) = Xall(:,:,k) - Xest;
                q = q + 1;
            end
        end
    end
    
    cvError(:,i) = sum(sum(errorTemp.^2,3));
end
figure;boxplot(cvError,Lambda);
xlabel('lambda');
ylabel('CV error');

%}



%{
sigma = .4;
cvError = zeros(n,length(Lambda));
for i = 1:length(Lambda)
    lambda = Lambda(i);
    param.lambda = lambda;
    alpha = mexLasso(X,D,param);
    Xest = D*alpha;
    
    Xpert1 = perturbImg(Y,width,height,'noise',sigma);
    Xtemp = Xpert1(ind,:);
    errorTemp1 = Xest - Xtemp;    
    
    Xpert2 = perturbImg(Y,width,height,'noise',sigma);
    Xtemp = Xpert2(ind,:);
    errorTemp2 = Xest - Xtemp;
   
    Xpert3 = perturbImg(Y,width,height,'noise',sigma);
    Xtemp = Xpert3(ind,:);
    errorTemp3 = Xest - Xtemp;

    
    Xpert4 = perturbImg(Y,width,height,'noise',sigma);
    Xtemp = Xpert4(ind,:);
    errorTemp4 = Xest - Xtemp;

    cvError(:,i) = (sum(errorTemp1.^2)+sum(errorTemp2.^2)+sum(errorTemp3.^2)+sum(errorTemp4.^2))';
end
figure;boxplot(cvError,Lambda);
xlabel('lambda');
ylabel('CV error');
%}