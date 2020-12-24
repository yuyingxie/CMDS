% Load relevant data:
%load('people2.mat'); Dis=people2; load('people2Labels.mat'); Labels = people2Labels';
%load('people3.mat'); Dis=people3; load('people3Labels.mat'); Labels = people3Labels';
load('people4.mat'); Dis=people4; load('people4Labels.mat'); Labels = people4Labels';
k=length(unique(Labels));

% Apply Double centering to distance matrix:
n = size(Dis,1);
J = eye(n) - ones(n,n)./n;
B = -J*(Dis.^2)*J;

% Compute eigendecomposition of B:
[V, D] = eig(B);
[D, I] = sort(diag(D), 'descend','ComparisonMethod','abs');
 V = V(:,I);
ExpRatio = D(1:end-1)./D(2:end);
eig_idx = find(D > 10^(-8));
r_ratio_est = find(max(ExpRatio(eig_idx(1:end-1)))==ExpRatio)
 
%% Create MDS embedding:
r = r_ratio_est;  %For People3 and 4
%r = 2; %For People2
Y = V(:,1:r)*sqrt(diag(D(1:r))); %CMDS Embedding
%  Z = linkage(Y,'single');
%  SL_clustering = cluster(Z,'Maxclust',k);
%  SLcounts=hist(SL_clustering,1:k)
kmeans_clustering = kmeans(Y, k,'Replicates',10);
kmeanscounts = hist(kmeans_clustering,1:k)
figure
gscatter(V(:,1)*sqrt(D(1)),V(:,2)*sqrt(D(2)),Labels,'bgmk','',[20 20 20 20],'off')
lgd = legend('Composers','Artists','Authors');
legend('Location','southeast');
lgd.FontSize = 16;
xlabel('First MDS Coordinate','fontsize',16)
ylabel('Second MDS Coordinate','fontsize',16)
%gscatter(Y(:,1),Y(:,2),Labels,'bgm','',[20 20 20],'off')
%gscatter(Y(:,1),Y(:,2),Labels,'bg','',[20 20],'off')
%axis square
%title('CMDS Embedding','fontsize',20)
 
  %% Quantify Accuracy
 M = confusionmat(Labels, kmeans_clustering);
 ClusterCounts = M(any(M,2),any(M,1))
 % Note: each row of ClusterCounts is a cancer, each column is a k-means
 % cluster
 
 ClassAccuracies = max(ClusterCounts,[],2)./sum(ClusterCounts,2)
 OverallAccuracy = sum(max(ClusterCounts,[],2))/sum(sum(ClusterCounts,2))
 OverallAccuracy_check = GetAccuracies(kmeans_clustering, Labels,k)
 
%  %% Compare with SLC
%  %Z = linkage(squareform(Dis),'single');
%  Z = linkage(squareform(Dis),'complete');
%  %Z = linkage(squareform(Dis),'average');
%  SL_clustering = cluster(Z,'MaxClust',k);
%  
%   M2 = confusionmat(Labels, SL_clustering);
%  ClusterCounts2 = M2(any(M2,2),any(M2,1))
%  % Note: each row of ClusterCounts is a cancer, each column is a k-means
%  % cluster
%  
%  ClassAccuracies2 = max(ClusterCounts2,[],2)./sum(ClusterCounts2,2)
%  OverallAccuracy2 = sum(max(ClusterCounts2,[],2))/sum(sum(ClusterCounts2,2))
 
 %% Approximate the SNR by calculating it in the full dimensional MDS embedding:
 
 Class_Idx{1,1} = 1:1:25;
 Class_Idx{1,2} = 26:1:50;
 Class_Idx{1,3} = 51:1:75;
 Class_Idx{1,4} = 76:1:100;
 
 pos_eig_idx = find(D>0);
 FullRankMDSEmbed = V(:,pos_eig_idx)*sqrt(diag(D(pos_eig_idx)));
 AllNoiseVecs = [];
 for i=1:k
     Mean{1,i} = mean(FullRankMDSEmbed(Class_Idx{1,i},:));
     NoiseVecs{1,i} = FullRankMDSEmbed(Class_Idx{1,i},:) - ones(25,1)*Mean{1,i};
     AllNoiseVecs = [AllNoiseVecs; NoiseVecs{1,i}]; 
 end
 MuDiffs = zeros(k,k);
for i=1:k
    for j=1:k
        MuDiffs(i,j) = norm(Mean{1,i}-Mean{1,j});
    end
end
[U,S,Vnoise] = svd(AllNoiseVecs/sqrt(n));
sigmaxsq = S(1,1)^2;
mudiff = MuDiffs(1,2); %this is actually minimal distance for all data sets, 2,3,4
SNR = mudiff^2/sigmaxsq

     
 
 
 
 
 