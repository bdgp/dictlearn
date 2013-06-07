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

numReplicates = 100;

sigma = 0;
numPatterns = 47;
Lambda = [0.01,0.05,0.1,0.2,0.3,0.5,0.8,1.0];
Gamma1 = [0.1,0.2,0.5,0.7,1,1.5,2,2.5,3];
initial_path = './32by16/randomStartParameters/';
for g = 1:length(Gamma1)
    gamma_char{g} = num2str(Gamma1(g));
end
for l = 1:length(Lambda)
    lambda_char{l} = num2str(Lambda(l));
end

doObj = 0
if doObj
    objFcnVal = zeros(length(Lambda),length(Gamma1),length(numPatterns));

    for k = 1:length(numPatterns)
        K = numPatterns(k);
        current_path0 = [initial_path, 'K=', num2str(K),'/']; 
        for l = 1:length(Lambda)
            
            lambda = Lambda(l);
            current_path1 = [current_path0, 'lambda=', num2str(lambda),'/']; 
            for g = 1:length(Gamma1)
                
                gamma1 = Gamma1(g);
                current_path2 = [current_path1, 'gamma1=', num2str(gamma1),'/']; 
                temp = dir([current_path2,'*.mat']);
                if length(temp) < 100
                    objFcnVal(l,g,k) = inf;
                else
                    load([current_path2,'bestDict.mat']);
                    objFcnVal(l,g,k) = min(R);
                end
            end
        end
    end

figure;imagesc(objFcnVal);colorbar;
set(gca,'YTick',1:length(lambda_char),'YTickLabel',lambda_char);
set(gca,'XTick',1:length(gamma_char),'XTickLabel',gamma_char);
xlabel('gamma value');
ylabel('lambda value');
title(['K = ',num2str(K),': objective function value']);
print(gcf,'-dpng',[current_path0,'objFcn.png']);


end
    
width = 32;
height = 16;
loadCV = 1;
if loadCV
    load('./32by16/data.mat');
    estStability = zeros(length(Lambda),length(Gamma1), ...
                         length(numPatterns));
    estDictError = estStability;

for k = 1:length(numPatterns)
    K = numPatterns(k);
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    for l = 1:length(Lambda)
        l = 8
        lambda = Lambda(l);

        current_path1 = [current_path0, 'lambda=', num2str(lambda),'/']; 
        for g = 1:length(Gamma1)
            g = 1
            gamma1 = Gamma1(g);
            current_path2 = [current_path1, 'gamma1=', num2str(gamma1),'/']; 
            temp = dir([current_path2,'*.mat']);
            if length(temp) < 100
                estStability(l,g,k) = inf;
                estDictError(l,g,k) = inf;
            else
                Dhat = zeros(size(X,1),K,numReplicates);

                for L = 1:numReplicates
                    load([current_path2,'rep',num2str(L), ...
                          'Dict.mat']);
                    Dhat(:,:,L) = D;
                end
                stop
                distMat = zeros(numReplicates,numReplicates);
                errorMat = zeros(numReplicates,numReplicates);
                for q = 1:numReplicates
                    for p = q:numReplicates
                        [distMat(q,p),distMat(p,q),errorMat(q,p),errorMat(p,q)] = distDictRegression(Dhat(:,:,q),Dhat(:,:,p));
                    end
                end
                figure;imagesc(distMat);
                colorbar;
                print(gcf,'-dpng',[current_path2,'distMatrixDict.png']);    
                save([current_path2,'distMatrixDict.mat'],'distMat','errorMat');    
                close all;
                estStability(l,g,k) = sum(distMat(:))/(numReplicates*(numReplicates-1));
                estDictError(l,g,k) = sum(errorMat(:))/(numReplicates*(numReplicates-1));
 
            end
            
        end
    end

    save([current_path0,'estStabDict.mat'],'numPatterns','estStability','estDictError');
    figure; imagesc(estStability); colorbar;
    set(gca,'YTick',1:length(lambda_char),'YTickLabel',lambda_char);
    set(gca,'XTick',1:length(gamma_char),'XTickLabel',gamma_char);
    xlabel('gamma value');
    ylabel('lambda value');
    title(['K = ',num2str(K),': estimation stability (uniqueness) of dictionaries']);
    print(gcf,'-dpng',[current_path0,'Dstability.png']);
   
    figure; imagesc(estDictError); colorbar;
    set(gca,'YTick',1:length(lambda_char),'YTickLabel',lambda_char);
    set(gca,'XTick',1:length(gamma_char),'XTickLabel',gamma_char);
    xlabel('gamma value');
    ylabel('lambda value');
    title(['K = ',num2str(K),': estimation stability (representation error) of dictionaries']);
    print(gcf,'-dpng',[current_path0,'Derror.png']);
   
    figure; scatter(estStability(:),estDictError(:));
    print(gcf,'-dpng',[current_path0,'stabVSError.png']);

end


end
%close all;

dictErrorAna = 0;
if dictErrorAna
    load('./32by16/data.mat');
    dictError = zeros(length(Lambda),length(Gamma1),length(numPatterns));
for k = 1:length(numPatterns)
    K = numPatterns(k);
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    for l = 1:length(Lambda)
        
        lambda = Lambda(l);
        current_path1 = [current_path0, 'lambda=', num2str(lambda),'/']; 
        for g = 1:length(Gamma1)
            
            gamma1 = Gamma1(g);
            current_path2 = [current_path1, 'gamma1=', num2str(gamma1),'/']; 
            temp = dir([current_path2,'*.mat']);
            if length(temp) < 100
                dictError(l,g,k) = inf;
            else
                Dhat = zeros(size(X,1),K,numReplicates);
                               
                for L = 1:numReplicates
                    load([current_path2,'rep',num2str(L), ...
                          'Dict.mat']);
                    Dhat(:,:,L) = D;
                end
                normConst = sum(Dhat(:).^2)/((numReplicates)*K);
                load([current_path2,'distMatrixDict.mat']);
                dictError(l,g,k) = sum(errorMat(:))/(numReplicates*(numReplicates-1)*normConst);
            end
        end
    end
    save([current_path0,'estStabDictError.mat'],'Lambda','Gamma1','dictError');

    figure; imagesc(dictError); colorbar;
    set(gca,'YTick',1:length(lambda_char),'YTickLabel',lambda_char);
    set(gca,'XTick',1:length(gamma_char),'XTickLabel',gamma_char);
    xlabel('gamma value');
    ylabel('lambda value');
    title(['K = ',num2str(K),': representation error of dictionaries']);
    print(gcf,'-dpng',[current_path0,'Derror.png']);
end
end