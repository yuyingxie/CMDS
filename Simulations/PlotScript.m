%% save emp_fail_prob_k2_dfixed_nscale.mat emp_fail_prob (use N)
%load('SuccessGrid_k2_d2_mudiff_small_svdHighRes.mat');
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(:,1:15)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(:,1:15),SuccessGrid(:,1:15),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',20)
ylabel('log(sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(1.951, 2.440, 20);
y2 = -33.48-x2;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS)) max(max(LogLogTotalSS)) min(min(LogSigmaSqGrid(:,1:15))) max(max(LogSigmaSqGrid(:,1:15)))])
colorbar

% %% save emp_fail_prob_k2_dfixed_nscale.mat emp_fail_prob (use n)
% load('SuccessGrid_k2_d2_mudiff_small_svdHighRes_BackupCopy.mat');
% figure
% %h=surf(LogLogNGrid(:,1:15),LogSigmaSqGrid(:,1:15),SuccessGrid(:,1:15),'FaceColor','inter')
% h=surf(log(log(2*exp(exp(LogLogNGrid(:,1:15))))),LogSigmaSqGrid(:,1:15),SuccessGrid(:,1:15),'FaceColor','inter')
% view(2)
% colorbar
% xlabel('log(log(n))','fontsize',20)
% ylabel('log(sigma^2)','fontsize',20)
% title('Empirical Probability of Exact Recovery','fontsize',16)
% 
% hold on
% x2 = linspace(1.951, 2.3680, 20);
% y2 = -33.55-x2;
% z = get(h,'ZData');
% set(h,'ZData',z-10)  
% plot(x2,y2,'m','linewidth',4)
% axis([min(min(LogLogNGrid(:,1:15))) max(max(LogLogNGrid(:,1:15))) min(min(LogSigmaSqGrid(:,1:15))) max(max(LogSigmaSqGrid(:,1:15)))])


%% 
load('SuccessGrid_k4_d2_mudiff_small_svd.mat');
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:17,1:20)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(1:17,1:20),SuccessGrid(1:17,1:20),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',20)
ylabel('log(sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(1.951, 2.440, 20);
y2 = -33.75-x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:17,1:20))) max(max(LogLogTotalSS(1:17,1:20))) min(min(LogSigmaSqGrid(1:17,1:20))) max(max(LogSigmaSqGrid(1:17,1:20)))])
colorbar


%% Plot for Fixed n, variable d
%load('SuccessGrid_k2_n100_ShiftSigma.mat');
figure
h=surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

%%
hold on
x=linspace(4.8,13,20);
y=-.5*x+.7;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSigmaSqGrid)) max(max(LogSigmaSqGrid))])