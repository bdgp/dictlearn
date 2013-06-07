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
numPatterns = 1:100;
initial_path = './32by16/randomStart/';

if ~exist('Load')
load('./32by16/data.mat');
ind = find(template(:,:,1)==1);
Load = 1;
end

rng('default');
[m,n] = size(X);
sigma = 0.1;
E = sigma*randn(m,n);
XNoise = X + E; % should leanr a dict from this...
XNoise(XNoise<0)=0;

stop
initial_path = './32by16/randomStart/';
rss = zeros(length(numPatterns),1);
for k = 1:length(numPatterns)
    K = numPatterns(k);
    current_path0 = [initial_path, 'K=', num2str(K),'/'];     
    load([current_path0,'bestDict.mat']);
    D = Dbest;
    lambda = 0;
    gamma1 = 0;

    
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

    
    alpha = mexLasso(XNoise,D,param);
    rss(k) = mean(sum((XNoise-D*alpha).^2))/m;
end

figure;plot(numPatterns,rss)
hline(sigma^2)
