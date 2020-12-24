N = 100;
d = 2.^(8:1:19);
CenteredNorm = zeros(1,length(d));
UncenteredNorm = zeros(1,length(d));

for i = 1:length(d)
    %Gaussian dis:
    %X = randn(N,d(i));
    %Uniform distribution
    X = sqrt(12)*rand(N,d(i))-sqrt(12)/2; 
    % mixture of Gaussians:
%     %%
%     U=rand(1,N);
%     U(find(U>.5))=3;
%     U(find(U<=.5))=-3;
%     X = (diag(U)*ones(N,d(i))+randn(N,d(i)))./sqrt(10);
    %%
    CenteredNorm(i) = norm(X*X'-d(i)*eye(N));
    UncenteredNorm(i) = norm(X*X');
end

%%
figure
scatter(log(d), log(CenteredNorm)) 

b = polyfit(log(d), log(CenteredNorm),1);
hold on
plot(log(d), b(1).*log(d)+b(2))
title('Centered Norm', 'fontsize',16)

slope_centered = b(1)

figure
scatter(log(d), log(UncenteredNorm)) 

b2 = polyfit(log(d), log(UncenteredNorm),1);
hold on
plot(log(d), b2(1).*log(d)+b2(2))
title('Uncentered Norm', 'fontsize',16)

slope_uncentered = b2(1)
    