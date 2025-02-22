% Set parameters (n and d, Covariance)
p=20;
n=floor(exp(exp(linspace(1.9,log(log(21000)),p))));
%n=floor(exp(exp(1.9)));
d = 20;  %Regime 1: fixed d
K=4; %Number of Nearest neighbord for KNN Cov
%Cov = MakeKNNCov(20,K);
load('KNNCov_d20_K4.mat');
SqrtCov = Cov^(1/2);
normCov = norm(Cov);
m=20;
sigmasq = exp(linspace(-79.8/normCov,-80.3/normCov,m));
%sigmasq = exp(-79.8/normCov);
sigma = sqrt(sigmasq);
J = 50; %Number of simulations
SuccessGrid = zeros(p,m);
LogSigmaSqGrid = zeros(p,m);
LogLogNGrid = zeros(p,m);


% Select Means:
M = [.0000001*eye(4) zeros(4,16)];
k = size(M,1);
r=k-1;


for s=1:p
    for q = 1:m
        N = k*n(s);
        gamma = d/N;
        Failure = ones(1,J); %AVL: Failure = zeros(1,J);
        for i=1:J
            X = zeros(N,d);
            for j=1:k
                X((j-1)*n(s)+1:j*n(s),:) = ones(n(s),1)*M(j,:)+sigma(q)*randn(n(s),d)*SqrtCov;
            end
            centeredX = X - ones(N,1)*mean(X);
            [V,S,R] = svd(centeredX);
            D = S.^2;
            %[V, D] = eig(centeredX*centeredX');
            [D, I] = sort(diag(D), 'descend');
             V = V(:,I);
             Z = linkage([real(V(:,1:r)*sqrt(diag(D(1:r)))) imag(V(:,1:r)*sqrt(diag(D(1:r))))],'single');
             SL_clustering = cluster(Z,'Maxclust',k);
             Failure(i) = 0;
             for j=1:k
                 if length(unique(SL_clustering((j-1)*n(s)+1:j*n(s))))>1
                     Failure(i)=1;
                 end
             end
% % AVL: when r=1, just use sign of top eigenvector for now as SL becomes comp bottleneck
%              posidx = find(V(:,1)>0);
%              negidx = find(V(:,1)<0);
%              if (isequal(negidx,(1:n(s))')||isequal(posidx,(1:n(s))'))==1
%                  Failure(i)=0;
%              end 
        end
        SuccessGrid(s,q) = 1-length(find(Failure==1))/J; %prob of exact recovery
        LogSigmaSqGrid(s,q) = log(normCov*sigmasq(q)); %Recods largest eigenvalue of the cov matrix
        LogLogNGrid(s,q) = log(log(n(s)));
        save('SuccessGrid_k4_d20_KNNCov.mat')
    end
end
% %% save emp_fail_prob_k2_dfixed_nscale.mat emp_fail_prob
% figure
% surf(LogLogNGrid,LogSigmaSqGrid,SuccessGrid)
% view(2)
% colorbar
% xlabel('log(log(n))')
% ylabel('log(sigma^2)')
% title('Probability of Exact Recovery')
% 
% save('SuccessGrid_k4_d10_GenCovHighRes.mat')