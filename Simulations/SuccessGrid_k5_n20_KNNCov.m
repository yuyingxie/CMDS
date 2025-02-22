% Set parameters (n and d)
n = 20;
p = 20; %number of d values (change to 8)
d=floor(exp(linspace(log(100),log(21000),p)));
m = 20; %number of sigma values
sigmasq = exp(linspace(-5.4,-6.2,m));
sigma = sqrt(sigmasq);
J = 20; %Number of simulations
SuccessGrid = zeros(p,m);
LogSigmaSqGrid = zeros(p,m);
LogDGrid = zeros(p,m);
NormCovGrid = zeros(p,m);
TraceCovGrid = zeros(p,m);
LogEffDimGrid = zeros(p,m);
K = 20; %Number of nearest neighbors to use when making covariance


for s=1:p
    Means = [.5 0 0 0 0 zeros(1,d(s)-5); 0 .5 0 0 0 zeros(1,d(s)-5); 0 0 .5 0 0 zeros(1,d(s)-5); 0 0 0 .5 0 zeros(1,d(s)-5); 0 0 0 0 .5 zeros(1,d(s)-5)];
    k = size(Means,1);
    r=k-1;
    M = [Means]; %augment with zeros
    N = k*n;
    for q = 1:m
        Failure = zeros(1,J);
        for i=1:J
            X = zeros(N,d(s));
            Cov = MakeKNNCov(d(s),K);
            SqrtCov = Cov^(1/2);
            NormCovGrid(s,q) = norm(Cov);
            TraceCovGrid(s,q) = trace(Cov);
            LogEffDimGrid(s,q) = log(TraceCovGrid(s,q)/NormCovGrid(s,q));
            for j=1:k
                X((j-1)*n+1:j*n,:) = ones(n,1)*M(j,:)+sigma(q)*randn(n,d(s))*SqrtCov;
            end
            centeredX = X - ones(N,1)*mean(X);
            [V, D] = eig(centeredX*centeredX');
            [D, I] = sort(diag(D), 'descend');
             V = V(:,I);
             Z = linkage([real(V(:,1:r))*sqrt(diag(D(1:r))) imag(V(:,1:r))*sqrt(diag(D(1:r)))],'single');
             SL_clustering = cluster(Z,'Maxclust',k);
             Failure(i) = 0;
             for j=1:k
                 if length(unique(SL_clustering((j-1)*n+1:j*n)))>1
                     Failure(i)=1;
                 end
             end
        end
        SuccessGrid(s,q) = 1-length(find(Failure==1))/J; %prob of exact recovery
        LogSigmaSqGrid(s,q) = log(sigmasq(q)*NormCovGrid(s,q));
        LogDGrid(s,q) = log(d(s));
        save('SuccessGrid_k5_n20_KNNCov.mat')
    end
end
% %% save emp_fail_prob_k2_dfixed_nscale.mat emp_fail_prob
% figure
% surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
% view(2)
% colorbar
% xlabel('log(d)')
% ylabel('log(sigma^2)')
% title('Probability of Exact Recovery')

%save('SuccessGrid_k5_n20_KNNCov.mat')