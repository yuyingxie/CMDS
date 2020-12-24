% n=20; %points per cluster
% sigma=.3;
% X1 = ones(n,1)*[10 0]+sigma*randn(n,2);
% X2 = ones(n,1)*[-4 .5]+sigma*randn(n,2);
% X3 = ones(n,1)*[-6 -.5]+sigma*randn(n,2);
% X=[X1; X2; X3];
load('AffectOfEigPertData.mat')
labels = [ones(n,1); 2*ones(n,1); 3*ones(n,1)];
figure
scatter(X(:,1),X(:,2),[],labels,'filled')
title('Original Data')
axis equal

%% Weight the eigenvalues equally:
centeredX = X - ones(3*n,1)*mean(X);
[V, D] = eig(centeredX*centeredX');
[D, I] = sort(diag(D), 'descend');
V = V(:,I);
figure
scatter(sqrt(D(1))*V(:,1),sqrt(D(2))*V(:,2),[],labels,'filled')
axis equal
title('MDS Embedding')
figure
scatter(V(:,1),V(:,2),[],labels,'filled')
axis equal
title('MDS Embedding without Weighting Eigenvalues')

%% Run k-means on true MDS embedding
MDS_kmeans = kmeans([sqrt(D(1))*V(:,1),sqrt(D(2))*V(:,2)],3);
UnweightedEigsMDS_kmeans = kmeans([V(:,1),V(:,2)],3);
figure
scatter(sqrt(D(1))*V(:,1),sqrt(D(2))*V(:,2),[],MDS_kmeans,'filled')
axis equal
title('K-means Classification MDS Embedding')
figure
scatter(V(:,1),V(:,2),[],UnweightedEigsMDS_kmeans,'filled')
axis equal
title('K-means Classification MDS Embedding without Weighting Eigenvalues')
