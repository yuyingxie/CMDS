function [C] = MakeKNNCov(d,K,s)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Generate points from 2-d uniform distribution in the cube:
X = s*rand(d,2);
[idx, D] = knnsearch(X,X,'K',K+1);
C = zeros(d,d);
for i=1:d
    for j=1:K+1
        if idx(i,j)==i
            C(i, idx(i,j))=1;
        else
            C(i, idx(i,j)) = D(i,j);
            C(idx(i,j), i) = D(i,j);
        end
    end
end
end

