d = 100;
N = 2.^(7:1:13);
CenteredNorm = zeros(1,length(N));
UncenteredNorm = zeros(1,length(N));

for i = 1:length(N)
    %Gaussian dis:
    X = randn(N(i),d);
    %Uniform distribution
    %X = sqrt(12)*rand(N(i),d)-sqrt(12)/2; 
    % mixture of Gaussians:
%     %%
%     U=rand(1,N);
%     U(find(U>.5))=3;
%     U(find(U<=.5))=-3;
%     X = (diag(U)*ones(N,d(i))+randn(N,d(i)))./sqrt(10);
    %%
    CenteredNorm(i) = norm(X*X'-d*eye(N(i)));
    UncenteredNorm(i) = norm(X*X');
end

%%
figure
scatter(log(N), log(CenteredNorm)) 

b = polyfit(log(N), log(CenteredNorm),1);
hold on
plot(log(N), b(1).*log(N)+b(2))
title('Centered Norm', 'fontsize',16)

slope_centered = b(1)

figure
scatter(log(N), log(UncenteredNorm)) 

b2 = polyfit(log(N), log(UncenteredNorm),1);
hold on
plot(log(N), b2(1).*log(N)+b2(2))
title('Uncentered Norm', 'fontsize',16)

slope_uncentered = b2(1)
    