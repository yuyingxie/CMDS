p=5;
M = 200000; %number of simulations
N = 10; %dimension of g_s
Sum_X_power_p = zeros(N,N);
Mean_X_power_p = zeros(N,N);
for i=1:M
    g = randn(N,1);
    Sum_X_power_p = Sum_X_power_p + (g*g'-eye(N))^p;
end
Mean_X_power_p = Sum_X_power_p/M;
figure
imagesc(Mean_X_power_p)
colorbar