%% Sim 1(aa)
load('SuccessGrid_k2_d2_mudiff_small_svdHighRes.mat');
mudiffsq = .0000002^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(:,1:15)))));
h=surf(LogLogTotalSS,LogSNRGrid(:,1:15),SuccessGrid(:,1:15),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%colormap gray
%title('Empirical Probability of Exact Recovery','fontsize',25)

hold on
x2 = linspace(1.951, 2.440, 20);
y2 = 2.63+x2;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS)) max(max(LogLogTotalSS)) min(min(LogSNRGrid(:,1:15))) max(max(LogSNRGrid(:,1:15)))])
axis square
colorbar


%% Sim 1(a)
load('SuccessGrid_k4_d2_mudiff_small_svd_SigmaZoom.mat');
mudiffsq = (sqrt(2)*.0000001)^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:17,1:20)))));
h=surf(LogLogTotalSS,LogSNRGrid(1:17,1:20),SuccessGrid(1:17,1:20),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(1.951, 2.440, 20);
y2 = 2.15+x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:17,1:20))) max(max(LogLogTotalSS(1:17,1:20))) min(min(LogSNRGrid(1:17,1:20))) max(max(LogSNRGrid(1:17,1:20)))])
axis square
colorbar

%% Sim. 1(b) - Toeplitz 

load('SuccessGrid_k4_d10_GenCovHighRes.mat');
mudiffsq = (sqrt(2)*.0000001)^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
row_idx = 17;
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:row_idx,:)))));
h=surf(LogLogTotalSS,LogSNRGrid(1:row_idx,:),SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(2, 2.36, 20);
y2 = .22+x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:row_idx,:))) max(max(LogLogTotalSS(1:row_idx,:))) min(min(LogSNRGrid(1:row_idx,:))) max(max(LogSNRGrid(1:row_idx,:)))])
axis square
colorbar

%% Sim 1(c) - KNN fixed d
load('SuccessGrid_k4_d20_KNNCov.mat');
mudiffsq = (sqrt(2)*.0000001)^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
row_idx = 17;
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:row_idx,:)))));
h=surf(LogLogTotalSS,LogSNRGrid(1:row_idx,:),SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x2 = linspace(2, 2.36, 20);
y2 = 1.44+x2;
z_max = max(max(get(h,'Zdata'))); 
%set(h,'ZData',z-10)  
line(x2,y2,z_max*ones(1,length(x2)),'linewidth',4,'Color','m')
axis([min(min(LogLogTotalSS(1:row_idx,:))) max(max(LogLogTotalSS(1:row_idx,:))) min(min(LogSNRGrid(1:row_idx,:))) max(max(LogSNRGrid(1:row_idx,:)))])
axis square
colorbar


%% Plot for Fixed n, variable d
% Sim 2(a)
load('SuccessGrid_k2_n100_ShiftSigma.mat');
mudiffsq = (sqrt(2))^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid,LogSNRGrid,SuccessGrid,'FaceColor','inter')
view(2)
axis square
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x=linspace(4.8,15.7,20);
y=.5*x-0;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSNRGrid)) max(max(LogSNRGrid))])

%% Sim 2(b)
load('SuccessGrid_k5_n20_rhoOne.mat');
mudiffsq = (sqrt(2)*.5)^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid,LogSNRGrid,SuccessGrid,'FaceColor','inter')
view(2)
axis square
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x=linspace(4.8,15.7,20);
y=.5*x+1;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSNRGrid)) max(max(LogSNRGrid))])

%% Sim 2(c)
load('SuccessGrid_k5_n20_GenCovHighRes.mat');
row_idx = 17;
mudiffsq = (sqrt(2)*.5)^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;

figure
h2=surf(LogDGrid(1:17,:),LogSNRGrid(1:17,:),SuccessGrid(1:17,:),'FaceColor','inter')
view(2)
axis square
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

x=linspace(4.5,8.9,20);
y=.5*x-.1;
%y=-.42*x-2;
z_max = max(max(get(h2,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid(1:row_idx,:))) max(max(LogDGrid(1:row_idx,:))) min(min(LogSNRGrid(1:row_idx,:))) max(max(LogSNRGrid(1:row_idx,:)))])

%% Sim 2(d) 

load('SuccessGrid_k5_n20_KNNCovNormalized_K10.mat');
row_idx=19;
mudiffsq = (sqrt(2)*.5)^2;
%LogSigmaSqGrid = ones(20,1)*logsigmasq; %Try removing norm cov, as it's too random
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid(1:row_idx,:),LogSNRGrid(1:row_idx,:), SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
axis square
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x=linspace(4.4,10,20);
y=.5*x+1;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid(1:row_idx,:))) max(max(LogDGrid(1:row_idx,:))) min(min(LogSNRGrid(1:row_idx,:))) max(max(LogSNRGrid(1:row_idx,:)))])



%% Prelim Sim 2(e) (k=3)
load('SuccessGrid_k3_n20_rhoLargeHighRes.mat');
mudiffsq = (sqrt(.4^2+.6^2))^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid,LogSNRGrid,SuccessGrid,'FaceColor','inter')
view(2)
axis square
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x1=linspace(4.6,9.5,20);
y1=.5*x1+1.2;
x2=linspace(9.5,16,20);
%y2=.75*x2-1.1750; %best fit
y2 =.78*x2-1.46; %to get same slope for 2(c) and 2(d)
x=[x1 x2];
y=[y1 y2];
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSNRGrid)) max(max(LogSNRGrid))])

%% Sim 2(f) (k=5)
load('SuccessGrid_k5_n20_rhoLarge.mat');
mudiffsq = (sqrt(.49^2+.51^2))^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid,LogSNRGrid,SuccessGrid,'FaceColor','inter')
view(2)
axis square
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)
hold on
x1=linspace(4.6,9,20);
y1=.5*x1+1.7;
x2=linspace(9,15,20);
%y2=.8*x2-1;  %best fit
y2=.78*x2-0.82;
x=[x1 x2];
y=[y1 y2];
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSNRGrid)) max(max(LogSNRGrid))])


%% Simulation 2(g): 2(e) but with Denoise Eigs
load('SuccessGrid_k3_n20_rhoLarge_DenoiseEigs.mat');
mudiffsq = (sqrt(.4^2+.6^2))^2;
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid,LogSNRGrid,SuccessGrid,'FaceColor','inter')
view(2)
axis square
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)

x=linspace(4.6,16,20);
y=.5*x+.5;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid)) max(max(LogDGrid)) min(min(LogSNRGrid)) max(max(LogSNRGrid))])


%% Sim 4(a) - Prelim

%load('SuccessGrid_k4_d10_BandedCov.mat');
row_idx=10;
figure
LogLogTotalSS = log(log(2*exp(exp(LogLogNGrid(1:row_idx,:)))));
h=surf(LogLogTotalSS,LogSigmaSqGrid(1:row_idx,:),SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
xlabel('log(log(N))','fontsize',30)
ylabel('log(\sigma_{max}^2)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

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
xlabel('log(d)','fontsize',30)
ylabel('log(\sigma_{max}^2)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

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




%% Sim 5(b) - fixed n
load('SuccessGrid_k5_n20_KNNCov.mat');
row_idx=19;
mudiffsq = (sqrt(2)*.5)^2;
%LogSigmaSqGrid = ones(20,1)*logsigmasq; %Try removing norm cov, as it's too random
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid(1:row_idx,:),LogSNRGrid(1:row_idx,:), SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)


% LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
% figure
% h=surf(LogDGrid(1:row_idx,:),LogSNRGrid(1:row_idx,:),SuccessGrid(1:row_idx,:),'FaceColor','inter')
% view(2)
% colorbar
% xlabel('log(d)','fontsize',20)
% ylabel('log(SNR)','fontsize',20)
% title('Empirical Probability of Exact Recovery','fontsize',16)

% hold on
% x=linspace(4.4,10,20);
% y=.4*x+2;
% z_max = max(max(get(h,'Zdata')));
% %set(h,'ZData',z-10)  
% line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
% axis([min(min(LogDGrid(1:row_idx,:))) max(max(LogDGrid(1:row_idx,:))) min(min(LogSNRGrid(1:row_idx,:))) max(max(LogSNRGrid(1:row_idx,:)))])


%% Sim 5(c) - fixed n
load('SuccessGrid_k5_n20_KNNCovNormalized.mat');
row_idx=19;
mudiffsq = (sqrt(2)*.5)^2;
%LogSigmaSqGrid = ones(20,1)*logsigmasq; %Try removing norm cov, as it's too random
LogSNRGrid = log(mudiffsq)-LogSigmaSqGrid;
figure
h=surf(LogDGrid(1:row_idx,:),LogSNRGrid(1:row_idx,:), SuccessGrid(1:row_idx,:),'FaceColor','inter')
view(2)
colorbar
xlabel('log(d)','fontsize',30)
ylabel('log(SNR)','fontsize',30)
%title('Empirical Probability of Exact Recovery','fontsize',16)

hold on
x=linspace(4.4,10,20);
y=.6*x-.2;
z_max = max(max(get(h,'Zdata')));
%set(h,'ZData',z-10)  
line(x,y,z_max*ones(1,length(x)),'linewidth',4,'Color','m')
axis([min(min(LogDGrid(1:row_idx,:))) max(max(LogDGrid(1:row_idx,:))) min(min(LogSNRGrid(1:row_idx,:))) max(max(LogSNRGrid(1:row_idx,:)))])
