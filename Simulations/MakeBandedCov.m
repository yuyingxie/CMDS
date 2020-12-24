function [Cov] = MakeBandedCov(d,NN,a)
%Input:
% d: dimension
% NN: number of neighbors to have dependecies with on each side (2*NN off diag bands) 
% a: the off-diagonal weight for the nonzero elements
Cov = sparse(d,d);
Cov = Cov + sparse(eye(d));
for i=1:NN
    Cov = Cov + sparse(diag(a*ones(1,d-i),-i))+sparse(diag(a*ones(1,d-i),i));
end
end

