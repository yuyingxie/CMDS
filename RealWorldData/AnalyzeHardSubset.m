% Analyze 3 kidney cancer types:   12    13    14
load('PANCANDigitLabels.mat'); Labels = PANCANDigitLabels;
load('PANCANFeatures.mat')

% Choose 3 types of kidney cancer:
index12 = find(PANCANDigitLabels==12); Labels(index12)=1;
index13 = find(PANCANDigitLabels==13); Labels(index13)=2;
index14 = find(PANCANDigitLabels==14); Labels(index14)=3;
SubsetIndex = [index12; index13; index14];

X = PANCANFeatures(SubsetIndex,:);
Labels = Labels(SubsetIndex);

% %% Denoise the data (otherwise we will just end up with outlier clusters in the embedding)
% K = 10;
% [idx, kNND] = knnsearch(X,X,'K',10);
% 
% %% Remove points with large kNN:
% figure
% scatter(1:size(X,1), sort(kNND(:,K)),'filled')
% cutoff = 400; %by inspection of graph
% indexToKeep = find(kNND(:,K)<cutoff);
% denoisedX = X(indexToKeep,:);
% denoisedLabels = Labels(indexToKeep);

%% Which genes to keep:

GenesToKeep = 1:20530; % just keep all of them

%% Just keep all points:

indexToKeep=1:size(X,1);
denoisedX = X(indexToKeep,GenesToKeep);
denoisedLabels = Labels(indexToKeep);

%% Calculate MDS matrix

%Center the Data:
meanX = mean(denoisedX);
N = size(denoisedX,1);
centeredX = denoisedX - ones(N,1)*mean(denoisedX);
%centeredX = centeredX./sqrt(Var(GenesToKeep)); %Optional: normalize the columns, all genes var 1
[V, D] = eig(centeredX*centeredX');
[D, I] = sort(diag(D), 'descend');
V = V(:,I);
ExpRatio = D(1:end-1)./D(2:end);
eig_idx = find(D > 10^(-8));
r_ratio_est = find(max(ExpRatio(eig_idx(1:end-1)))==ExpRatio)
 
 %% Implement SL Clustering
 r=r_ratio_est; %Calculate using eigenratio huer: r=3.
 k = 3; 
 Y = V(:,1:r)*sqrt(diag(D(1:r))); %CMDS Embedding
 Z = linkage(Y,'single');
 SL_clustering = cluster(Z,'Maxclust',k);
 %SLcounts=hist(SL_clustering,1:k)
 kmeans_clustering = kmeans(Y, k,'Replicates',20);
 kmeanscounts = hist(kmeans_clustering,1:k);
 
 %% Quantify Accuracy
 M = confusionmat(denoisedLabels, kmeans_clustering);
 ClusterCounts = M(any(M,2),any(M,1));
 % Note: each row of ClusterCounts is a cancer, each column is a k-means
 % cluster
 
 %ClassAccuracies = max(ClusterCounts,[],2)./sum(ClusterCounts,2)
 %OverallAccuracy = sum(max(ClusterCounts,[],2))/sum(sum(ClusterCounts,2))
 
 OverallAccuracy_kmeans = GetAccuracies(kmeans_clustering, denoisedLabels,k)
 OverallAccuracy_SL = GetAccuracies(SL_clustering, denoisedLabels,k)
 
 %% Plot CMDS Embedding
figure
%scatter(Y(:,1),Y(:,2),[],Labels,'filled')
gscatter(Y(:,1),Y(:,2),Labels,'bgmk','',[20 20 20 20],'off')
axis square
% lgd = legend('Breast invasive carcinoma','Kidney renal cc carcinoma','Lung adenocarcinoma', 'Thyroid carcinoma')
% legend('Location','southeast')
% lgd.FontSize = 14;
% xlabel('First MDS Coordinate','fontsize',16)
% ylabel('Second MDS Coordinate','fontsize',16)
% xlim([-200, 200])
 %xlabel('First multidimensional scaling coordinate','fontsize',16)
 %ylabel('Second multidimensional scaling coordinate','fontsize',16)
 %title('CMDS Embedding of TCGA PAN Cancer Data','fontsize',16)
 
%  figure
%  scatter(Y(:,3),Y(:,4),[],Labels,'filled')
%  
%  figure
%  scatter(Y(:,5),Y(:,6),[],Labels,'filled')

%% Calc the SNR for this example
mu = zeros(4,length(GenesToKeep));
mu(1,:) = mean(PANCANFeatures(index12,GenesToKeep));
mu(2,:) = mean(PANCANFeatures(index13,GenesToKeep)); 
mu(3,:) = mean(PANCANFeatures(index14,GenesToKeep));  

MuDiffs = zeros(3,3);
for i=1:3
    for j=1:3
        MuDiffs(i,j) = norm(mu(i,:)-mu(j,:));
    end
end

%%
mudiff = MuDiffs(2,3);

CenteredX1 = PANCANFeatures(index12,GenesToKeep) - ones(length(index12),1)*(mu(1,:));
CenteredX2 = PANCANFeatures(index13,GenesToKeep) - ones(length(index13),1)*(mu(2,:));
CenteredX3 = PANCANFeatures(index14,GenesToKeep) - ones(length(index14),1)*(mu(3,:));

% % For calculating sigma_max, info reported below:
NoiseVecs = [CenteredX1; CenteredX2; CenteredX3];
[U,S,Vnoise] = svd(NoiseVecs/sqrt(N));
sigmaxsq = S(1,1)^2;

%sigmaxsq = 3.3975e+03;
SNR = mudiff^2/sigmaxsq

% %% Run k-means on full data set
%  kmeans_clusteringAllData = kmeans(centeredX, k,'Replicates',10);
%  
% %Quantify Accuracy of k-means on full data
%  M2 = confusionmat(denoisedLabels, kmeans_clusteringAllData);
%  ClusterCounts2 = M2(any(M2,2),any(M2,1))
%  % Note: each row of ClusterCounts is a cancer, each column is a k-means
%  % cluster
%  
%  ClassAccuracies2 = max(ClusterCounts2,[],2)./sum(ClusterCounts2,2)
%  OverallAccuracy2 = sum(max(ClusterCounts2,[],2))/sum(sum(ClusterCounts2,2))
 
