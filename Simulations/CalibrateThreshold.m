FirstEigenvalue = zeros(1000,1);
SecondEigenvalues = zeros(1000,1);
Ratio = zeros(1000,1);
N = 1000;
d = 200;
for i=1:1000
    W = randn(N,d);
    centeredW = W - ones(N,1)*mean(W);
    [V, D] = eig(centeredW'*centeredW);
    [D, I] = sort(diag(D), 'descend');
    V = V(:,I);
    Ratio(i) = D(1)/D(2);
end
prctile(Ratio,98)