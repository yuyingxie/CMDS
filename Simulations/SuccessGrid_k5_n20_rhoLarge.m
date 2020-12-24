% Set parameters (n and d)
n = 20;
dmin = 100;
p = 20; %number of d values (change to 8)
%d=floor(exp(linspace(log(100),log(700000),p)));
d=floor(exp(linspace(log(100),log(6553600),p)));
%sigmasqmin = .0003;
m = 30; %number of sigma values
%sigmasq = exp(linspace(log(.003),log(.1928),m));
sigmasq = exp(linspace(-11,-4,m));
sigma = sqrt(sigmasq);
%n = [60 120 240 480 960];
%sigmasq=.01:.07:.4;
%p = size(n,2); % number of sample sizes
%m = size(sigma,2); % number of sigma values
J = 20; %Number of simulations
SuccessGrid = zeros(p,m);
LogSigmaSqGrid = zeros(p,m);
LogDGrid = zeros(p,m);

% Select Means (first 6 coordinates):
% Data set 1 (k=2, colinear):

% Data set 2 (k=6, r=5):
% Means = eye(6);
% Data set 3 (k=6, r=1):
% Means = [-3 0 0 0 0 0; -2 0 0 0 0 0; -1 0 0 0 0 0; 1 0 0 0 0 0; 2 0 0 0 0 0; 3 0 0 0 0 0];

for s=1:p
    Means = [0 0 zeros(1,d(s)-2); .49 .51 zeros(1,d(s)-2); -.51 -.49 zeros(1,d(s)-2); 0 0 .49 .51 zeros(1,d(s)-4); 0 0 -.51 -.49 zeros(1,d(s)-4)];
    k = size(Means,1);
    r=k-1;
    M = [Means]; %augment with zeros
    N = k*n;
    for q = 1:m
        Failure = zeros(1,J);
        for i=1:J
            X = zeros(N,d(s));
            for j=1:k
                X((j-1)*n+1:j*n,:) = ones(n,1)*M(j,:)+sigma(q)*randn(n,d(s));
            end
            centeredX = X - ones(N,1)*mean(X);
            [V, D] = eig(centeredX*centeredX');
            [D, I] = sort(diag(D), 'descend');
             V = V(:,I);
             Z = linkage(V(:,1:r)*sqrt(diag(D(1:r))),'single');
             SL_clustering = cluster(Z,'Maxclust',k);
             Failure(i) = 0;
             for j=1:k
                 if length(unique(SL_clustering((j-1)*n+1:j*n)))>1
                     Failure(i)=1;
                 end
             end
        end
        SuccessGrid(s,q) = 1-length(find(Failure==1))/J; %prob of exact recovery
        LogSigmaSqGrid(s,q) = log(sigmasq(q));
        LogDGrid(s,q) = log(d(s));
        save('SuccessGrid_k5_n20_rhoLarge.mat')
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

%save('SuccessGrid_k3_n20_rhoLargeHighRes.mat')