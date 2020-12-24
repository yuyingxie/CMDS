%% Sim 1(a)
load('SuccessGrid_k2_d2_mudiff_small_svdHighRes.mat');
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(:,1:15)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(:,1:15),SuccessGrid(:,1:15),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',20)
ylabel('log(\sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(1.951, 2.440, 20);
y2 = -33.48-x2;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS)) max(max(LogLogTotalSS)) min(min(LogSigmaSqGrid(:,1:15))) max(max(LogSigmaSqGrid(:,1:15)))])
colorbar


%% Sim 1(b)
load('SuccessGrid_k4_d2_mudiff_small_svd_SigmaZoom.mat');
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:17,1:20)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(1:17,1:20),SuccessGrid(1:17,1:20),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',20)
ylabel('log(\sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(1.951, 2.440, 20);
y2 = -33.7-x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:17,1:20))) max(max(LogLogTotalSS(1:17,1:20))) min(min(LogSigmaSqGrid(1:17,1:20))) max(max(LogSigmaSqGrid(1:17,1:20)))])
colorbar


%% Plot for Fixed n, variable d
% Sim 2(a)
load('SuccessGrid_k2_n100_ShiftSigma.mat');
figure
h=surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(\sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x=linspace(4.8,15.7,20);
y=-.5*x+.7;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSigmaSqGrid)) max(max(LogSigmaSqGrid))])

%% Sim 2(b)
load('SuccessGrid_k5_n20_rhoOne.mat');
figure
h=surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(\sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x=linspace(4.8,15.7,20);
y=-.5*x-1.7;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSigmaSqGrid)) max(max(LogSigmaSqGrid))])

%% Prelim Sim 2(c) (k=3)
load('SuccessGrid_k3_n20_rhoLargeHighRes.mat');
figure
h=surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(\sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x1=linspace(4.8,9.5,20);
y1=-.5*x1-1.8;
x2=linspace(9.5,15,20);
y2=-.8*x2+1.05;
x=[x1 x2];
y=[y1 y2];
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSigmaSqGrid)) max(max(LogSigmaSqGrid))])

%% Sim 2(d) (k=5)
load('SuccessGrid_k5_n20_rhoLarge.mat');
figure
h=surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(\sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x1=linspace(4.6,9,20);
y1=-.5*x1-2.4;
x2=linspace(9,15,20);
y2=-.8*x2+.3;
x=[x1 x2];
y=[y1 y2];
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSigmaSqGrid)) max(max(LogSigmaSqGrid))])

%% Sim 3(a) - Prelim

load('SuccessGrid_k4_d10_GenCov.mat');
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:9,:)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(1:9,:),SuccessGrid(1:9,:),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',20)
ylabel('log(\sigma_{max}^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(2, 2.36, 20);
y2 = -31.825-x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:9,:))) max(max(LogLogTotalSS(1:9,:))) min(min(LogSigmaSqGrid(1:9,:))) max(max(LogSigmaSqGrid(1:9,:)))])
colorbar

%% Sim. 3(a) 

load('SuccessGrid_k4_d10_GenCovHighRes.mat');
row_idx = 17;
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:row_idx,:)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(1:row_idx,:),SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',20)
ylabel('log(\sigma_{max}^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(2, 2.36, 20);
y2 = -31.77-x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:row_idx,:))) max(max(LogLogTotalSS(1:row_idx,:))) min(min(LogSigmaSqGrid(1:row_idx,:))) max(max(LogSigmaSqGrid(1:row_idx,:)))])
colorbar

%% Sim 3(b)
%load('SuccessGrid_k5_n20_GenCovHighRes.mat');
row_idx = 17;

% % Effective Dimension:
% figure
% h=surf(LogEffDimGrid(1:row_idx,:),LogSigmaSqGrid(1:row_idx,:),SuccessGrid(1:row_idx,:),'FaceColor','inter')
% view(2)
% colorbar
% xlabel('log(Eff. Dim.)','fontsize',20)
% ylabel('log(\sigma_{max}^2)','fontsize',20)
% title('Empirical Probability of Exact Recovery','fontsize',16)
% 
% x=linspace(2.8,7.1,20);
% y=-.5*x-1.45;
% %y=-.42*x-2;
% z_max = max(max(get(h,'Zdata')));
% %set(h,'ZData',z-10)  
% line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
% axis([min(min(LogEffDimGrid(1:row_idx,:))) max(max(LogEffDimGrid(1:row_idx,:))) min(min(LogSigmaSqGrid(1:row_idx,:))) max(max(LogSigmaSqGrid(1:row_idx,:)))])

% Regular Dimension:
figure
h2=surf(LogDGrid(1:17,:),LogSigmaSqGrid(1:17,:),SuccessGrid(1:17,:),'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(\sigma_{max}^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

x=linspace(4.5,8.9,20);
y=-.5*x-.58;
%y=-.42*x-2;
z_max = max(max(get(h2,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid(1:row_idx,:))) max(max(LogDGrid(1:row_idx,:))) min(min(LogSigmaSqGrid(1:row_idx,:))) max(max(LogSigmaSqGrid(1:row_idx,:)))])

%% Sim 4(a) - Prelim

%load('SuccessGrid_k4_d10_BandedCov.mat');
row_idx=10;
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:row_idx,:)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(1:row_idx,:),SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',20)
ylabel('log(\sigma_{max}^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(2, 2.36, 20);
y2 = -31.825-x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:row_idx,:))) max(max(LogLogTotalSS(1:row_idx,:))) min(min(LogSigmaSqGrid(1:row_idx,:))) max(max(LogSigmaSqGrid(1:row_idx,:)))])
colorbar


%% Sim 4(b)
load('SuccessGrid_k5_n20_BandedCov.mat')

%Real Dim:
figure
h=surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(\sigma_{max}^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)

x=linspace(4.5,9.5,20);
y=-.36*x-1.82;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSigmaSqGrid)) max(max(LogSigmaSqGrid))])

% % Effective dim
% 
% figure
% h2=surf(LogEffDimGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
% view(2)
% colorbar
% xlabel('log(Eff Dim)','fontsize',20)
% ylabel('log(\sigma_{max}^2)','fontsize',20)
% title('Empirical Probability of Exact Recovery','fontsize',16)
% 
% x=linspace(2.5,9.5,20);
% y=-.36*x-2.45;
% z_max = max(max(get(h2,'Zdata')));
% %set(h,'ZData',z-10)  
% line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
% axis([min(min(LogEffDimGrid)) max(max(LogEffDimGrid)) min(min(LogSigmaSqGrid)) max(max(LogSigmaSqGrid))])

%% Sim 5(b)
load('SuccessGrid_k5_n20_KNNCov.mat');
figure
h=surf(LogDGrid,LogSigmaSqGrid,SuccessGrid,'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',20)
ylabel('log(\sigma^2)','fontsize',20)
title('Empirical Probability of Exact Recovery','fontsize',16)
