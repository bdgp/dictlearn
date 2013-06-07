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
    %template = generateTemplate();
    %template = imresize(template,0.5,'nearest');
    ind = find(template(:,:,1)==1);
    Load = 1;
end
%rng(215);
%indSet = testTrainIndProduce(size(X,2),5,'kmeans');
%holdInd = indSet{1};
numReplicates = 2;

numPatterns = 1:2:100;
%Lambda = [0,0.1,0.2,0.5,1];
%Gamma1 = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];

initial_path = ['./32by16/CV/',num2str(numReplicates),'fold/'];
loadCV = 1;
if loadCV
rss = zeros(length(numPatterns),1);
numNonZeros = rss;
estStability = rss;
for k = 1:length(numPatterns)
    K = numPatterns(k);    
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    
    rss(k) = 0;
    %Xhat = zeros(size(X,1),length(holdInd),numReplicates);
    for L = 1:numReplicates
        load([current_path0,'cvRep',num2str(L), ...
              '.mat']);
        rss(k) = rss(k) + testError;
                %alpha = mexLasso(X(:,testInd),DTrain,param);
                %alpha = mexLasso(X(:,holdInd),DTrain,param);
                %Xhat(:,:,L) = DTrain*alpha;

                %numNonZerosTemp = sum(sum(abs(alpha)>1e-4))/length(testInd); 
                %numNonZeros(k,l,g) = numNonZeros(k,l,g) + numNonZerosTemp;
                % should add model complexity as a criterion
                % or delete pixels in a nice way
                % or leave out one group of sample to do stab
    end
            
    %{
    estStabVarTemp = zeros(1,length(holdInd));
    Xbar = mean(Xhat,3);
    for q = 1:length(holdInd)
                
        for L = 1:numReplicates
            sqDistTemp = sum((Xhat(:,q,L) - Xbar(:,q)).^2);
            estStabVarTemp(q) = estStabVarTemp(q)+sqDistTemp;
        end
                
    end
    estStabVarTemp = estStabVarTemp/numReplicates;
    estStabMagTemp = sum(Xbar.^2);
            %estStabMagTemp = sum(X(:,holdInd).^2);
           
            estStabTemp = estStabVarTemp./estStabMagTemp;
            estStability(k,l,g) = mean(estStabTemp); % can also median
                                                     
            %rss(k,l,g) = rss(k,l,g)/L;
            %numNonZeros(k,l,g) = numNonZeros(k,l,g)/L;
            %}
end
end

figure; plot(numPatterns,rss);
xlabel('K');
ylabel('test error');
print(gcf,'-dpng',['32by16/CV/',num2str(numReplicates),'fold/rss.png']);

stop
for k = 1:5
    K = numPatterns(k);
    
    
    current_path0 = [initial_path, 'K=', num2str(K),'/']; 
    
    escvTemp = reshape(estStability(k,:,:),length(Lambda),length(Gamma1));
    figure;imagesc(escvTemp);colorbar;
    set(gca,'YTick',1:length(lambda_char),'YTickLabel',lambda_char);
    set(gca,'XTick',1:length(gamma_char),'XTickLabel',gamma_char);
    xlabel('gamma value');
    ylabel('lambda value');
    title(['K = ',num2str(K),': estimation stability cross-validation']);
    set(gcf,'PaperPositionMode','auto'); % preserve the size
    print(gcf,'-dpng',[current_path0,'ESCV_img.png']);
    
    %{
    rssTemp = reshape(rss(k,:,:),length(Lambda),length(Gamma1));
    figure;imagesc(log(rssTemp));colorbar;
    set(gca,'YTick',1:length(lambda_char),'YTickLabel',lambda_char);
    set(gca,'XTick',1:length(gamma_char),'XTickLabel',gamma_char);
    xlabel('gamma value');
    ylabel('lambda value');
    title(['K = ',num2str(K),': cross-validation residual sum of squares']);
    set(gcf,'PaperPositionMode','auto'); % preserve the size
    print(gcf,'-dpng',[current_path0,'CV_RSS.png']);
    
    
    
    numTemp = reshape(numNonZeros(k,:,:),length(Lambda),length(Gamma1));
    figure;imagesc(numTemp);colorbar;
    set(gca,'YTick',1:length(lambda_char),'YTickLabel',lambda_char);
    set(gca,'XTick',1:length(gamma_char),'XTickLabel',gamma_char);
    ylabel('lambda value');
    xlabel('gamma value');
    title(['K = ',num2str(K),': number of nonzeros']);
    set(gcf,'PaperPositionMode','auto'); % preserve the size
    print(gcf,'-dpng',[current_path0,'CV_sparsityAlphs.png']);
    %}
    
    %set(gcf,'PaperPositionMode','auto'); % preserve the size
    %set(gcf,'InvertHardcopy','off'); % preserve the background color
    %print(gcf,'-dpng',[current_path0,'CV_RSS.png']);
    %input('press enter to continue');
end
stop
%{
% find similar patterns across the CV partitions
for k = 1:3
  K = KK(k);  
  current_path0 = [initial_path, 'K=', num2str(K),'/']; 
      
  for l = 1:length(Lambda)
    lambda = Lambda(l);
    current_path1 = [current_path0, 'lambda=',num2str(lambda),'/']; 
      
    for g = 1:length(Gamma1)
      gamma1 = Gamma1(g);
      current_path2 = [current_path1, 'gamma1=', num2str(gamma1),'/'];
      current_path3 = [current_path2,'crossValidation/'];
      load([current_path2,'spca_Result.mat']);
      figure('Position',[1423,17,1244,901]);
      Dmain = D;
      for i = 1:10
          subplot(10,11,11*(i-1)+1);
          displayTI(Dmain(:,i),tf.t,tf.p_scale);
      end
      input('press enter to continue');
      for L = 1:10
          load([current_path3,'spca_cv_partition',num2str(L), ...
                '.mat']);
          Dtemp = D;
          for i = 1:10
              corrTemp = corr(Dmain(:,i),Dtemp);
              maxIndex = find(corrTemp == max(corrTemp));
              subplot(10,11,11*(i-1)+L+1);
              displayTI(Dtemp(:,maxIndex),tf.t,tf.p_scale);
              title(corrTemp(maxIndex));
          end
          input('press enter to continue');
      end
      stop
     end
   end
end
%}