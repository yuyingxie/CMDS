addpath(genpath('/Users/annalittle/MDS/Simulations'))

RunRNAMix = 1; %Choose 1 for RNA mix and 0 for cell mix

if RunRNAMix == 1
    % RNA mix (Linnorm):
    %load('RNAmix1.mat'); Data = RNAmix1'; Labels = RNAmix1Labels; mix1_idx=find(Labels==1); Labels(mix1_idx) = 2; Labels=Labels-1;
    load('RNAmix2.mat'); Data = RNAmix2'; Labels = RNAmix2Labels; mix1_idx=find(Labels==1); Labels(mix1_idx) = 2; Labels=Labels-1;
elseif RunRNAMix == 0
    % Cell mix (Linnorm):
    load('/Users/annalittle/Dropbox/MSU/scRNAseq/DataSets/sc10x5clLinnorm.mat');load('/Users/annalittle/Dropbox/MSU/scRNAseq/DataSets/truelabelssc10x5clnopreprocessing.mat'); Data=sc10x5clLinnorm'; Labels= truelabelssc10x5clnopreprocessing;
end

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

GT1_idx = find(Labels==1);
GT2_idx = find(Labels==2);
GT3_idx = find(Labels==3);
GT4_idx = find(Labels==4);
GT5_idx = find(Labels==5);
if RunRNAMix==1
    GT6_idx = find(Labels==6);
    GT7_idx = find(Labels==7);
end

%% Which genes to keep:

GenesToKeep = 1:size(Data,2); % just keep all of them

% sorted_var = sort(var(Data),'descend');
% cutoff = sorted_var(2000+1);
% GenesToKeep = find((var(Data)>cutoff));

X = Data(:,GenesToKeep);

%% Calculate MDS matrix

%Center the Data:
meanX = mean(X);
N = size(X,1);
centeredX = X - ones(N,1)*mean(X);
%centeredX = centeredX./sqrt(Var(GenesToKeep)); %Optional: normalize the columns, all genes var 1
[V, D] = eig(centeredX*centeredX');
[D, I] = sort(diag(D), 'descend');
V = V(:,I);
ExpRatio = D(1:end-1)./D(2:end);
eig_idx = find(D > 10^(-8));
%r_ratio_est = find(max(ExpRatio(eig_idx(1:end-1)))==ExpRatio)
r_ratio_est = 2
 
 %% Implement SL Clustering
 r=r_ratio_est; %Calculate using eigenratio huer: r=2.
 k = length(unique(Labels)); 
 Y = V(:,1:r)*sqrt(diag(D(1:r))); %CMDS Embedding
 Z = linkage(Y,'single');
 SL_clustering = cluster(Z,'Maxclust',k);
 %SLcounts=hist(SL_clustering,1:k)
 kmeans_clustering = kmeans(Y, k,'Replicates',20);
 kmeanscounts = hist(kmeans_clustering,1:k);
 
 %% Quantify Accuracy
 M = confusionmat(Labels, kmeans_clustering);
 ClusterCounts = M(any(M,2),any(M,1));
 % Note: each row of ClusterCounts is a cancer, each column is a k-means
 % cluster
 
 %ClassAccuracies = max(ClusterCounts,[],2)./sum(ClusterCounts,2)
 %OverallAccuracy = sum(max(ClusterCounts,[],2))/sum(sum(ClusterCounts,2))
 
 OverallAccuracy_kmeans = GetAccuracies(kmeans_clustering, Labels,k)
 OverallAccuracy_SL = GetAccuracies(SL_clustering, Labels,k)
 
 %% Plot CMDS Embedding
figure
if RunRNAMix == 1
    gscatter(Y(:,1),Y(:,2),Labels,'bgmkyrc','',[20 20 20 20],'off')
    figure
    gscatter(Y(:,1),Y(:,2),SL_clustering,'bgmkyrc','',[20 20 20 20],'off')
    title('SL Clustering')
elseif RunRNAMix == 0
    gscatter(Y(:,1),Y(:,2),Labels,'bgmky','',[20 20 20 20],'off')
end
    
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
mu = zeros(k,length(GenesToKeep));
mu(1,:) = mean(X(GT1_idx,:));
mu(2,:) = mean(X(GT2_idx,:));   
mu(3,:) = mean(X(GT3_idx,:)); 
mu(4,:) = mean(X(GT4_idx,:));
mu(5,:) = mean(X(GT5_idx,:)); 
if RunRNAMix==1
    mu(6,:) = mean(X(GT6_idx,:));
    mu(7,:) = mean(X(GT7_idx,:));
end

MuDiffs = zeros(k,k);
for i=1:k
    for j=1:k
        MuDiffs(i,j) = norm(mu(i,:)-mu(j,:));
    end
end

%%
mudiff = min(MuDiffs(MuDiffs>0));

CenteredX1 = X(GT1_idx,:) - ones(length(GT1_idx),1)*(mu(1,:));
CenteredX2 = X(GT2_idx,:) - ones(length(GT2_idx),1)*(mu(2,:));
CenteredX3 = X(GT3_idx,:) - ones(length(GT3_idx),1)*(mu(3,:));
CenteredX4 = X(GT4_idx,:) - ones(length(GT4_idx),1)*(mu(4,:));
CenteredX5 = X(GT5_idx,:) - ones(length(GT5_idx),1)*(mu(5,:));
if RunRNAMix==1
    CenteredX6 = X(GT6_idx,:) - ones(length(GT6_idx),1)*(mu(6,:));
    CenteredX7 = X(GT7_idx,:) - ones(length(GT7_idx),1)*(mu(7,:));
end

% % For calculating sigma_max, info reported below:
if RunRNAMix==1
    NoiseVecs = [CenteredX1; CenteredX2; CenteredX3; CenteredX4; CenteredX5; CenteredX6; CenteredX7];
elseif RunRNAMix==0
    NoiseVecs = [CenteredX1; CenteredX2; CenteredX3; CenteredX4; CenteredX5];
end
[U,S,Vnoise] = svd(NoiseVecs/sqrt(N));
sigmaxsq = S(1,1)^2;

SNR = mudiff^2/sigmaxsq

% %% Run k-means on full data set
%  kmeans_clusteringAllData = kmeans(centeredX, k,'Replicates',20);
%  
% %Quantify Accuracy of k-means on full data
%  M2 = confusionmat(Labels, kmeans_clusteringAllData);
%  ClusterCounts2 = M2(any(M2,2),any(M2,1))
%  
%  OverallAccuracy_fulldata = GetAccuracies(kmeans_clusteringAllData, Labels,k)
 
% %% Create umap embedding for comparison:
% addpath(genpath('umapFileExchange (1.5.2)'))
% 
% [reduction, umap, clusterIdentifiers, extras]=run_umap(X);
% 
% figure
% if RunRNAMix == 1
%     gscatter(reduction(:,1),reduction(:,2),Labels,'bgmkyrc','',[20 20 20 20],'off')
% elseif RunRNAMix == 0
%     gscatter(reduction(:,1),reduction(:,2),Labels,'bgmky','',[20 20 20 20],'off')
% end

%% Create tsne embedding for comparison:

tsne_coordinates = tsne(X);
figure
if RunRNAMix == 1
    gscatter(tsne_coordinates(:,1),tsne_coordinates(:,2),Labels,'bgmkyrc','',[20 20 20 20],'off')
elseif RunRNAMix == 0
    gscatter(tsne_coordinates(:,1),tsne_coordinates(:,2),Labels,'bgmky','',[20 20 20 20],'off')
end