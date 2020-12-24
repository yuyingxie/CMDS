% Record of parameters chosen and set-up
n = 200;
%d = 1000; %for paper example
d = 200; 
Means = [2 1; -1 -1; 1 -1; -2 1; 0 0];
k = size(Means,1);
N = k*n;
sigma = .25;
M = [Means zeros(k,d-2)]; %augment with zeros

X = [ones(n,1)*M(1,:)+sigma*randn(n,d); 
     ones(n,1)*M(2,:)+sigma*randn(n,d); 
     ones(n,1)*M(3,:)+sigma*randn(n,d);
     ones(n,1)*M(4,:)+sigma*randn(n,d);
     ones(n,1)*M(5,:)+sigma*randn(n,d)];

GT = [ones(1,n) 2*ones(1,n) 3*ones(1,n) 4*ones(1,n) 5*ones(1,n)];

centeredX = X - ones(N,1)*mean(X);
[V, D] = eig(centeredX*centeredX');
[D, I] = sort(diag(D), 'descend');
 V = V(:,I);
 % % SL clustering on MDS embedding:
 %Z = linkage(V(:,1:r)*sqrt(diag(D(1:r))),'single');
 %clustering = cluster(Z,'Maxclust',k);
 
 %% Estimate r:

ExpRatio = D(1:end-1)./D(2:end);
ratio_r_est = find(max(ExpRatio(1:100))==ExpRatio)
ExpDiff = D(1:end-1)-D(2:end);
maxgap_r_est = find(max(ExpDiff(1:100))==ExpDiff)
 
 %% Cluster using ratio
 
 r=200;
 %r=ratio_r_est;
 MDS_clustering = kmeans(V(:,1:r)*sqrt(diag(D(1:r))),k,'replicates',10);
 M = confusionmat(GT, MDS_clustering);
 ClusterCounts = M(any(M,2),any(M,1));
 % Note: each row of ClusterCounts is a cancer, each column is a k-means
 % cluster
 
 ClassAccuracies_ratio = max(ClusterCounts,[],2)./sum(ClusterCounts,2)
 OverallAccuracy_ratio = sum(max(ClusterCounts,[],2))/sum(sum(ClusterCounts,2))
 
  %% Cluster using max gap
 
 r=maxgap_r_est;
 MDS_clustering = kmeans(V(:,1:r)*sqrt(diag(D(1:r))),k,'replicates',10);
 M = confusionmat(GT, MDS_clustering);
 ClusterCounts = M(any(M,2),any(M,1));
 % Note: each row of ClusterCounts is a cancer, each column is a k-means
 % cluster
 
 ClassAccuracies_maxgap = max(ClusterCounts,[],2)./sum(ClusterCounts,2)
 OverallAccuracy_maxgap = sum(max(ClusterCounts,[],2))/sum(sum(ClusterCounts,2))


%%
figure
scatter(V(:,1),V(:,2),25,GT,'filled')
title('2D CMDS Embedding of the Data','fontsize',16)

% %% Compute Spectral Embedding of the data:
% Dis = zeros(N,N);
% for i=1:N
%     for j=1:N
%         Dis(i,j) = norm(X(i,:)-X(j,:));
%     end
% end
% 
% b=prctile(Dis(:),5);
% W = exp(-Dis.^2/b^2);
% 
% Deg = diag(sum(W,1));
% 
% SymLap = eye(N) - (Deg^(-.5))*W*(Deg^(-.5));
% 
% [V2, D2] = eig(SymLap);
% [D2, I2] = sort(diag(D2), 'ascend');
%  V2 = V2(:,I2);
%  
%  figure
% %scatter(V2(:,1),V2(:,2),25,GT,'filled')
% scatter(V2(:,2),V2(:,3),25,GT,'filled')
% title('Spectral Clustering Embedding','fontsize',16)

% %% Hierarchical Clustering directly on distances
% 
%  %% Compare with SLC
%  %Z = linkage(squareform(Dis),'single');
%  %Z = linkage(X,'single');
%  %Z = linkage(X,'average');
%  %Z = linkage(X,'centroid');
%  %Z = linkage(X,'median');
%  %Z = linkage(X,'weighted');
%  %Z = linkage(X,'complete');
%  Z = linkage(X,'ward');
%  %Z = linkage(squareform(Dis),'average');
%  SL_clustering = cluster(Z,'MaxClust',k);
%  
%  M2 = confusionmat(GT, SL_clustering);
%  ClusterCounts2 = M2(any(M2,2),any(M2,1))
%  % Note: each row of ClusterCounts is a cancer, each column is a k-means
%  % cluster
%  
%  ClassAccuracies2 = max(ClusterCounts2,[],2)./sum(ClusterCounts2,2)
%  OverallAccuracy2 = sum(max(ClusterCounts2,[],2))/sum(sum(ClusterCounts2,2))
%  
%  figure
%  scatter(X(:,1),X(:,2),25,SL_clustering,'filled')
% title('Hierarchical Clustering Results','fontsize',16)
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



 
    