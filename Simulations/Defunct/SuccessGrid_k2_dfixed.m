% Set parameters (n and d)
%n = [60 120 240 480 960 1920];
nmin = 12800;
p = 4; %number of n values (change to 8)
n=nmin*((2.^(0:p-1)));
sigmasqmin = .0003;
m = 6; %number of sigma values
sigmasq = sigmasqmin*((4.^(0:5)));
sigma = sqrt(sigmasq);
%n = [60 120 240 480 960];
%sigmasq=.01:.07:.4;
%p = size(n,2); % number of sample sizes
%m = size(sigma,2); % number of sigma values
J = 10; %Number of simulations
d = 2;  %Regime 1: fixed d
SuccessGrid = zeros(p,m);
LogSigmaSqGrid = zeros(p,m);
LogLogNGrid = zeros(p,m);

% Select Means (first 6 coordinates):
% Data set 1 (k=2, colinear):
Means = [1 0; 0 1];
% Data set 2 (k=6, r=5):
% Means = eye(6);
% Data set 3 (k=6, r=1):
% Means = [-3 0 0 0 0 0; -2 0 0 0 0 0; -1 0 0 0 0 0; 1 0 0 0 0 0; 2 0 0 0 0 0; 3 0 0 0 0 0];
k = size(Means,1);
M = [Means]; %augment with zeros

for s=1:p
    for q = 1:m
        N = k*n(s);
        gamma = d/N;
        Failure = zeros(1,J);
        for i=1:J
            X = zeros(N,d);
            for j=1:k
                X((j-1)*n(s)+1:j*n(s),:) = ones(n(s),1)*M(j,:)+sigma(q)*randn(n(s),d);
            end
            centeredX = X - ones(N,1)*mean(X);
            [V, D] = eig(centeredX*centeredX');
            [D, I] = sort(diag(D), 'descend');
             V = V(:,I);
             Z = linkage(V(:,1:k),'single');
             SL_clustering = cluster(Z,'Maxclust',k);
             Failure(i) = 0;
             for j=1:k
                 if length(unique(SL_clustering((j-1)*n(s)+1:j*n(s))))>1
                     Failure(i)=1;
                 end
             end
        end
        SuccessGrid(s,q) = 1-length(find(Failure==1))/J; %prob of exact recovery
        LogSigmaSqGrid(s,q) = log(sigmasq(q));
        LogLogNGrid(s,q) = log(log(n(s)));
        save('SuccessGrid_k2_d2_nto102400.mat')
    end
end
%% save emp_fail_prob_k2_dfixed_nscale.mat emp_fail_prob
figure
surf(LogLogNGrid,LogSigmaSqGrid,SuccessGrid)
view(2)
colorbar
xlabel('log(log(n))')
ylabel('log(sigma^2)')
title('Probability of Exact Recovery')

save('SuccessGrid_k2_d2_nto102400.mat')