% Let's try to establish the effective dimension for the KNN example we had
% before: K=20, a=1
K=20;
s=1;
d=100;
C=MakeKNNCov(d,K,s);
[V D] = eig(C);
D = sort(diag(D));
maxeig = max(abs(D));
effdim = sqrt(sum(abs(D).^4)/maxeig^4)