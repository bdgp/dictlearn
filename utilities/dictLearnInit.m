function D = dictLearnInit(X,K,method,p)
% X: the data matrix
% K: the number of dictinary elements
% method: which initialization method to use?
% p: (optional) how many elements in cluster are used to compute
% the initial values?

if ~exist('p')
    p = 0; % default to take the cluster centroids
end

if strmatch(method,'kmeans','exact')
    loadPath = ['~/fruitfly/image/osDict/principalPatternPaper/cluster/kmeans/K=',num2str(K),'/clusterResult.mat'];
    load(loadPath);
end

if strmatch(method,'nmf','exact')
    loadPath = ['~/fruitfly/image/osDict/principalPatternPaper/dictLearn/bestDict/randomStart/K=',num2str(K),'/bestDict.mat'];
    load(loadPath);
    template = generateTemplate();
    template = imresize(template,0.5,'nearest');
    D = Dbest(template(:,:,1)==1,:);
    return;
end

if strmatch(method,'random','exact')
    ind = randsample(1:size(X,2),K);
    D = X(:,ind);
    return;
end

D = zeros(size(X,1),K);
for i = 1:K
    indTemp = find(idx == i);
    clusterSize = length(indTemp);
    % no checking for zero cluster size!!
    if p == 0        
        selectInd = 1:clusterSize; % select everything inside the cluster
    else
        selectInd = randsample(clusterSize,min(p,clusterSize));
    end
    Xtemp = X(:,indTemp(selectInd));
    D(:,i) = mean(Xtemp,2);
end
    
        
        
