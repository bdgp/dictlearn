sigma = 0;
numPatterns = 1:2:100;
numReplicates = 100;
dictError = zeros(1,length(numPatterns));
initial_path = ['./32by16/randomStartSigma=',num2str(sigma),'/'];
initial_path = ['./32by16/randomStart/'];

for k = 1:length(numPatterns)
    K = numPatterns(k);
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    load([current_path0,'distMatrixDict.mat']);
    dictError(k) = sum(errorMat(:))/(numReplicates*(numReplicates-1))/K;
end

figure;plot(numPatterns,dictError);
title(['sigma=',num2str(sigma)]);