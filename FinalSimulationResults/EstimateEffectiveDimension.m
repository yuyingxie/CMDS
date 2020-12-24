% Let's try to establish the effective dimension for the KNN example we had
% before: K=20, a=1
K=20;
s=1;
d=2210;
C=MakeKNNCov(d,K,s);
[V D] = eig(C);
D = sort(diag(D));
maxeig = max(abs(D));
effdim = sum(abs(D).^2)/maxeig^2