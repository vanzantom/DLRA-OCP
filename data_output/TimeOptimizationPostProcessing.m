clear all;close all;
timeOtP=zeros(4,2);
timePtO=zeros(4,2);
load('TimeOptimization_Nh_225_N_128.mat');
timefull(1)=Full.time;
timeOtP(:,1)=OtP.time;
timePtO(:,1)=PtO.time;
nvec(1)=128;mvec(1)=225;
sizeproblem(1)=mvec(1)*nvec(1);
load('TimeOptimization_Nh_961_N_256.mat');
timefull(2)=Full.time;
timeOtP(:,2)=OtP.time;
timePtO(:,2)=PtO.time;
nvec(2)=256;mvec(2)=961;
sizeproblem(2)=mvec(2)*nvec(2);
load('TimeOptimization_Nh_3969_N_512.mat');
timefull(3)=Full.time;
timeOtP(:,3)=OtP.time;
timePtO(:,3)=PtO.time;
nvec(3)=3969;mvec(3)=512;
sizeproblem(3)=mvec(3)*nvec(3);

figure(1)
loglog(sizeproblem,timefull,'k-d','Linewidth',2,'Markersize',10);
hold on
loglog(sizeproblem,timeOtP(1,:),'k-o','Linewidth',2,'Markersize',10);
loglog(sizeproblem,timeOtP(2,:),'k-x','Linewidth',2,'Markersize',10);
loglog(sizeproblem,timeOtP(3,:),'k-s','Linewidth',2,'Markersize',10);
loglog(sizeproblem,timeOtP(4,:),'k-*','Linewidth',2,'Markersize',10);

loglog(sizeproblem,timePtO(1,:),'k--o','Linewidth',2,'Markersize',10, 'HandleVisibility', 'off');
loglog(sizeproblem,timePtO(2,:),'k--x','Linewidth',2,'Markersize',10, 'HandleVisibility', 'off');
loglog(sizeproblem,timePtO(3,:),'k--s','Linewidth',2,'Markersize',10, 'HandleVisibility', 'off');
loglog(sizeproblem,timePtO(4,:),'k--*','Linewidth',2,'Markersize',10, 'HandleVisibility', 'off');

loglog(sizeproblem,10^(-1.5)*sizeproblem,'k','Linewidth',2)
loglog(sizeproblem,(mvec+nvec),'k-.','Linewidth',2)
grid on;
xlabel('Problem size (m\times n)')
ylabel('Computational time');
set(gca,'fontsize',16)
legend('Full-rank','r=4','r=8','r=16','r=32','m\times n','m+n','location','northwest')
xlim([1e4,5.1e6])
ylim([7e1,1e5])
yticks([1e2, 1e3, 1e4, 1e5])
filename='../fig_output/TimeOptimization.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 