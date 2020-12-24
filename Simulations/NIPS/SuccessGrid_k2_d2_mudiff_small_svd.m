% Set parameters (n and d)
%n = [60 120 240 480 960 1920 3840];
%n=floor(50*(1.41.^(1.41.^(0:8))));
n=floor(exp(exp(linspace(1.9,log(log(42300)),20))));
p=size(n,2);
%n = [60 120 240 480 960 1920];
%nmin = 12800;
%p = 4; %number of n values (change to 8)
%n=nmin*((2.^(0:p-1)));
%sigmasqmin = .00000000000000000005;
%sigmasqmin = .000000000000000005;
%m = 6; %number of sigma values
%sigmasq = sigmasqmin*((4.^(0:5)));
sigmasq = exp(linspace(-35.5,-36.5,20));
sigma = sqrt(sigmasq);
m=size(sigma,2);
%n = [60 120 240 480 960];
%sigmasq=.01:.07:.4;
%p = size(n,2); % number of sample sizes
%m = size(sigma,2); % number of sigma values
J = 20; %Number of simulations
d = 2;  %Regime 1: fixed d
SuccessGrid = zeros(p,m);
LogSigmaSqGrid = zeros(p,m);
LogLogNGrid = zeros(p,m);

% Select Means (first 6 coordinates):
% Data set 1 (k=2, colinear):
Means = [.0000001 0; 0 .0000001];
% Data set 2 (k=6, r=5):
% Means = eye(6);
% Data set 3 (k=6, r=1):
% Means = [-3 0 0 0 0 0; -2 0 0 0 0 0; -1 0 0 0 0 0; 1 0 0 0 0 0; 2 0 0 0 0 0; 3 0 0 0 0 0];
k = size(Means,1);
r=k-1;
M = [Means]; %augment with zeros

for s=1:p
    for q = 1:m
        N = k*n(s);
        gamma = d/N;
        Failure = ones(1,J); %AVL: Failure = zeros(1,J);
        for i=1:J
            X = zeros(N,d);
            for j=1:k
                X((j-1)*n(s)+1:j*n(s),:) = ones(n(s),1)*M(j,:)+sigma(q)*randn(n(s),d);
            end
            centeredX = X - ones(N,1)*mean(X);
            [V,S,R] = svd(centeredX);
            D = S.^2;
            %[V, D] = eig(centeredX*centeredX');
            [D, I] = sort(diag(D), 'descend');
             V = V(:,I);
% AVL: Just look at eigenvector for now as SL becomes comp bottleneck
%              Z = linkage(V(:,1:r)*diag(D(1:r)),'single');
%              SL_clustering = cluster(Z,'Maxclust',k);
%              Failure(i) = 0;
%              for j=1:k
%                  if length(unique(SL_clustering((j-1)*n(s)+1:j*n(s))))>1
%                      Failure(i)=1;
%                  end
%              end
             posidx = find(V(:,1)>0);
             negidx = find(V(:,1)<0);
             if (isequal(negidx,(1:n(s))')||isequal(posidx,(1:n(s))'))==1
                 Failure(i)=0;
             end 
        end
        SuccessGrid(s,q) = 1-length(find(Failure==1))/J; %prob of exact recovery
        LogSigmaSqGrid(s,q) = log(sigmasq(q));
        LogLogNGrid(s,q) = log(log(n(s)));
        save('SuccessGrid_k2_d2_mudiff_small_svdHighRes.mat')
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

save('SuccessGrid_k2_d2_mudiff_small_svdHighRes.mat')