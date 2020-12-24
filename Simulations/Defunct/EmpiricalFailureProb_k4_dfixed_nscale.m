% Set parameters (n and d)
%n = [30 60 120 240 480 960];
%n = [30 60 120];
n = [60 120 240];
J = 2000; %Number of simulations
d = 4;  %Regime 1: fixed d
% d = log(n); %Regime 2
% d = n; %Regime 3: fixed ratio
% d = n^2; %Regime 4: high-dimensional

% Select Means (first 6 coordinates):

% Data set 1 (k=2, colinear):
Means = eye(4);
% Data set 2 (k=6, r=5):
% Means = eye(6);
% Data set 3 (k=6, r=1):
% Means = [-3 0 0 0 0 0; -2 0 0 0 0 0; -1 0 0 0 0 0; 1 0 0 0 0 0; 2 0 0 0 0 0; 3 0 0 0 0 0];
k = size(Means,1);
M = [Means]; %augment with zeros

emp_fail_prob = zeros(1, size(n,1));
for s=1:size(n,2)
    N = k*n(s);
    gamma = d/N;
    if d < N
        sigma = 1/(20*max(sqrt(log(n(s))),sqrt(d)));
    else
        sigma = 1/(20*max(sqrt(gamma),d^(1/4)));
    end
    Failure = zeros(1,J);
    for i=1:J
        X = zeros(N,d);
        for j=1:k
            X((j-1)*n(s)+1:j*n(s),:) = ones(n(s),1)*M(j,:)+sigma*randn(n(s),d);
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
    emp_fail_prob(s) = length(find(Failure==1))/J;
end
save emp_fail_prob_k4_dfixed_nscale.mat emp_fail_prob