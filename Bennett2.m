p=10;
M = 1; %number of simulations
N1 = 10; %dimension of g_s
N2 = 15;
Sum_X_power_p = zeros(N1+N2,N1+N2);
Mean_X_power_p = zeros(N1+N2,N1+N2);
for i=1:M
    g1 = randn(N1,1);
    g2 = randn(N2,1);
    X = [zeros(N1,N1) g1*g2'; g2*g1' zeros(N2,N2)];
    Sum_X_power_p = Sum_X_power_p + (X)^p;
end
Mean_X_power_p = Sum_X_power_p/M;
figure
imagesc(Mean_X_power_p)
colorbar