% Analyze four most numerous cancer types:   5    13    16    30
load('PANCANDigitLabels.mat')
load('PANCANFeatures.mat')

index4 = find(PANCANDigitLabels==4);
index5 = find(PANCANDigitLabels==5);
index11 = find(PANCANDigitLabels==11);
index13 = find(PANCANDigitLabels==13);
index16 = find(PANCANDigitLabels==16);
index17 = find(PANCANDigitLabels==17);
index23 = find(PANCANDigitLabels==23);
index30 = find(PANCANDigitLabels==30);
index32 = find(PANCANDigitLabels==32);
SubsetIndex = [index4; index5; index11; index13; index16; index17; index23; index30; index32];

X = PANCANFeatures(SubsetIndex,:);
Labels = PANCANDigitLabels(SubsetIndex);

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

%% Only keep high variance genes:
Var = zeros(1, size(X,2));
Var5 = zeros(1, size(X,2));
Var13 = zeros(1, size(X,2));
Var16 = zeros(1, size(X,2));
Var30 = zeros(1, size(X,2));
for i=1:size(X,2)
    Var(i) = var(X(:,i));
    Var5(i) = var(PANCANFeatures(index5,i));
    Var13(i) = var(PANCANFeatures(index13,i));
    Var16(i) = var(PANCANFeatures(index16,i));
    Var30(i) = var(PANCANFeatures(index30,i));
end
VarClusterRatio = zeros(size(Var));
for i=1:length(Var)
    smallest_var = min([Var5(i) Var13(i) Var16(i) Var30(i)]);
    largest_var = max([Var5(i) Var13(i) Var16(i) Var30(i)]);
    VarClusterRatio(i) = largest_var/smallest_var;
end
% Keep genes with 1) high variance and 2) each cluster variance not too
% unbalanced.
%GenesToKeep = find(Var>12 & VarClusterRatio < 5); 99.5% acc on both
%MDS+cluster and cluster
%GenesToKeep = find(Var>7.15 & VarClusterRatio < 1.5); %25 genes; 87% acc MDS, 89% acc cluster
%GenesToKeep = find(Var>4.25 & VarClusterRatio < 1.5); %100 genes; 96.7% MDS, 97.9% full data

%GenesToKeep = find(Var>26); %this will give top 10 genes; 0.9664 accuracy with r=7
%GenesToKeep = find(Var>15.6); %this gives top 100 genes; 0.9963 accuracy with r=7
%GenesToKeep = find(Var>6.4806); % this will give top 1000 genes; 0.9946
%acc w/ r=7

%GenesToKeep = find(Var>0); %Keep all
GenesToKeep = find(Var>0 & VarClusterRatio < 3); %Keep all with same ratio;

%% Just keep all points:
indexToKeep=1:size(X,1);
denoisedX = X(indexToKeep,GenesToKeep);
denoisedLabels = Labels(indexToKeep);

%% Calculate MDS matrix

%Center the Data:
meanX = mean(denoisedX);
N = size(denoisedX,1);
centeredX = denoisedX - ones(N,1)*mean(denoisedX);
centeredX = centeredX./sqrt(Var(GenesToKeep)); %Optional: normalize the columns, all genes var 1
[V, D] = eig(centeredX*centeredX');
[D, I] = sort(diag(D), 'descend');
 V = V(:,I);
 
 %% Implement SL Clustering
 %r = 7; %From inspection of spectrum
 r=20;
 k = 9; 
 Y = V(:,1:r)*sqrt(diag(D(1:r))); %CMDS Embedding
 Z = linkage(Y,'single');
 SL_clustering = cluster(Z,'Maxclust',k);
 SLcounts=hist(SL_clustering,1:k)
 kmeans_clustering = kmeans(Y, k,'Replicates',20);
 kmeanscounts = hist(kmeans_clustering,1:k)
 
 %% Quantify Accuracy
 M = confusionmat(denoisedLabels, kmeans_clustering);
 ClusterCounts = M(any(M,2),any(M,1))
 % Note: each row of ClusterCounts is a cancer, each column is a k-means
 % cluster
 
 ClassAccuracies = max(ClusterCounts,[],2)./sum(ClusterCounts,2)
 OverallAccuracy = sum(max(ClusterCounts,[],2))/sum(sum(ClusterCounts,2))
 
 %% Plot CMDS Embedding
 figure
 scatter(Y(:,1),Y(:,2),[],Labels,'filled')
 xlabel('First CMDS Coordinate','fontsize',16)
 ylabel('Second CMDS Coordinate','fontsize',16)
 title('CMDS Embedding of TCGA PAN Cancer Data','fontsize',16)
 
%  figure
%  scatter(Y(:,3),Y(:,4),[],Labels,'filled')
%  
%  figure
%  scatter(Y(:,5),Y(:,6),[],Labels,'filled')

%% Calc the SNR for this example
mu = zeros(4,length(GenesToKeep));
mu(1,:) = mean(PANCANFeatures(index5,GenesToKeep));
mu(2,:) = mean(PANCANFeatures(index13,GenesToKeep)); 
mu(3,:) = mean(PANCANFeatures(index16,GenesToKeep));  
mu(4,:) = mean(PANCANFeatures(index30,GenesToKeep)); 

MuDiffs = zeros(4,4);
for i=1:4
    for j=1:4
        MuDiffs(i,j) = norm(mu(i,:)-mu(j,:));
    end
end

mudiff = MuDiffs(1,3);

CenteredX1 = PANCANFeatures(index5,GenesToKeep) - ones(length(index5),1)*(mu(1,:));
CenteredX2 = PANCANFeatures(index13,GenesToKeep) - ones(length(index13),1)*(mu(2,:));
CenteredX3 = PANCANFeatures(index16,GenesToKeep) - ones(length(index16),1)*(mu(3,:));
CenteredX4 = PANCANFeatures(index30,GenesToKeep) - ones(length(index30),1)*(mu(4,:));

NoiseVecs = [CenteredX1; CenteredX2; CenteredX3; CenteredX4];
sigmaxsq = norm(NoiseVecs'*NoiseVecs/size(NoiseVecs,1))

%% Run k-means on full data set
 kmeans_clusteringAllData = kmeans(centeredX, k,'Replicates',10);
 
  %% Quantify Accuracy of k-means on full data
 M2 = confusionmat(denoisedLabels, kmeans_clusteringAllData);
 ClusterCounts2 = M2(any(M2,2),any(M2,1))
 % Note: each row of ClusterCounts is a cancer, each column is a k-means
 % cluster
 
 ClassAccuracies2 = max(ClusterCounts2,[],2)./sum(ClusterCounts2,2)
 OverallAccuracy2 = sum(max(ClusterCounts2,[],2))/sum(sum(ClusterCounts2,2))