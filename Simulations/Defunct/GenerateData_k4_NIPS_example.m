% Record of parameters chosen and set-up
n = 200;
d = 100; 
Means = [1 1; -1 -1; 1 -1; -1 1; 0 0];
k = size(Means,1);
r = 2;
N = k*n;
sigma = .1;
M = [Means zeros(k,d-2)]; %augment with zeros

load('Data_d100_n200_k4_r2.mat');

centeredX = X - ones(N,1)*mean(X);
[V, D] = eig(centeredX*centeredX');
[D, I] = sort(diag(D), 'descend');
 V = V(:,I);
 Z = linkage(V(:,1:r),'single');
 SL_clustering = cluster(Z,'Maxclust',k);
 Failure = 0;
 for j=1:k
     if length(unique(SL_clustering((j-1)*n+1:j*n)))>1
         Failure=1;
     end
 end
Failure

%%
scatter(X(:,33),X(:,50),25,SL_clustering,'filled')
%xlabel('Data Dimension 33')
%ylabel('Data Dimension 50')
title('Coordinates 33 and 50 of Original Data','fontsize',16)

%%
figure
scatter(V(:,1),V(:,2),25,SL_clustering,'filled')
title('CMDS Embedding of the Data','fontsize',16)

 
    