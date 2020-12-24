% Record of parameters chosen and set-up
n = 200;
%d = 1000; %for paper example
d = 1000; 
Means = [1 1; -1 -1; 1 -1; -1 1; 0 0];
%Means = [2 1; -1 -1; 1 -1; -2 1; 0 0];
k = size(Means,1);
r = 2;
N = k*n;
%sigma = .3;
M = [Means zeros(k,d-2)]; %augment with zeros

% X = [ones(n,1)*M(1,:)+sigma*randn(n,d); 
%      ones(n,1)*M(2,:)+sigma*randn(n,d); 
%      ones(n,1)*M(3,:)+sigma*randn(n,d);
%      ones(n,1)*M(4,:)+sigma*randn(n,d);
%      ones(n,1)*M(5,:)+sigma*randn(n,d)];

%load('k5_n200_sig.2.mat');
%load('k5_n200_sig.25.mat')
load('k5_n200_sig.3.mat');
GT = [ones(1,n) 2*ones(1,n) 3*ones(1,n) 4*ones(1,n) 5*ones(1,n)];

centeredX = X - ones(N,1)*mean(X);
[V, D] = eig(centeredX*centeredX');
[D, I] = sort(diag(D), 'descend');
 V = V(:,I);
 Y = V(:,1:r)*sqrt(diag(D(1:r))); %CMDS embedding
 % % SL clustering on MDS embedding:
 %Z = linkage(V(:,1:r)*sqrt(diag(D(1:r))),'single');
 %clustering = cluster(Z,'Maxclust',k);
 MDS_clustering = kmeans(Y,k,'replicates',10);
 MDSAccuracy = GetAccuracies(MDS_clustering,GT',k)

%%
figure
scatter(X(:,33),X(:,50),25,GT,'filled')
axis square
xlim([1.05*min(X(:,33)) 1.05*max(X(:,33))])
ylim([1.05*min(X(:,50)) 1.05*max(X(:,50))])
%xlabel('Data Dimension 33')
%ylabel('Data Dimension 50')
%title('Coordinates 33 and 50 of Original Data','fontsize',16)

%%
figure
scatter(Y(:,1),Y(:,2),25,GT,'filled')
axis square
%title('CMDS Embedding of the Data','fontsize',16)

%% Hierarchical Clustering directly on distances

 %% Compare with SLC
 %Z = linkage(squareform(Dis),'single');
 %Z = linkage(X,'single');
 %Z = linkage(X,'average');
 %Z = linkage(X,'centroid');
 %Z = linkage(X,'median');
 %Z = linkage(X,'weighted');
 %Z = linkage(X,'complete');
 Z = linkage(X,'ward');
 %Z = linkage(squareform(Dis),'average');
 SL_clustering = cluster(Z,'MaxClust',k);
 HierAccuracy = GetAccuracies(SL_clustering,GT',k)
 
 
 figure
 scatter(X(:,1),X(:,2),25,SL_clustering,'filled')
 axis square
%title('Hierarchical Clustering Results','fontsize',16)

% Plot the accuracy of MDS + k-means as a function of r
% (note: for each r values, 50 simulations were run and accuracy averaged)
load('VaryRSim_Sigma3.mat')
figure
errorbar(R, ClusterAccuracy, ClusterAccuracySE,'LineWidth',2)
ylim([.75 1])
xlim([2 500])
axis square
box off

% % 
% %  figure
% %  %dendrogram(Z,'Labels',GT)
% % %  [H,T] = dendrogram(Z,'colorthreshold',29.14);
% 
% %% For fun let's try PCA:
% 
% [V3, D3] = eig(centeredX'*centeredX);
% [D3, I3] = sort(diag(D3), 'descend');
%  V3 = V3(:,I3);
%  
%  pc1 = centeredX*V3(:,1);
%  pc2 = centeredX*V3(:,2);
%  
% figure
% scatter(pc1,pc2,25,GT,'filled')
% title('PCA Embedding of the Data','fontsize',16)

% %% Examine ratio of Eigenvalues:
% 
% ExpRatio = D(1:end-1)./D(2:end);
% ratio_r_est = find(max(ExpRatio(1:100))==ExpRatio)
% ExpDiff = D(1:end-1)-D(2:end);
% maxgap_r_est = find(max(ExpDiff(1:100))==ExpDiff)
% 

 
    