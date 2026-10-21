clear all;close all;
load('EffectiveCurvature_Nh_961_N_500_Deltat_1.00e-03.mat');

figure(1)
step=100;
semilogy(tvec,kappaeff(:,1),'k-o','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
hold on
semilogy(tvec,1./sigmar(:,1),'k--o','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,kappaeff(:,2),'k-x','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,1./sigmar(:,2),'k--x','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,kappaeff(:,3),'k-s','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,1./sigmar(:,3),'k--s','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,kappaeff(:,4),'k-*','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,1./sigmar(:,4),'k--*','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
grid on;
legend('r=4','','r=8','','r=16','','r=32','');xlabel('t')
title('$$\kappa_{eff}(t)\quad and\quad 1/\sigma_r(Y(t))$$','interpreter','latex');
set(gca,'fontsize',18)
filename='../fig_output/effective_curvature.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 

figure(2)
semilogy(tvec,Fperp(:,1),'k-o','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
hold on
semilogy(tvec,Fperp(:,2),'k-x','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,Fperp(:,3),'k-s','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,Fperp(:,4),'k-*','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
grid on;
legend('r=4','r=8','r=16','r=32');
xlabel('t');ylabel('$$\|\mathcal{P}_Y^\perp F(t,Y,z)\|_{M,\Xi}$$','interpreter','latex')
set(gca,'fontsize',18)
yticks([1e-5 1e-3])
filename='../fig_output/residualF.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 

figure(3)
semilogy(tvec,FxQperp(:,1),'k-o','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
hold on
semilogy(tvec,FxQperp(:,2),'k-x','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,FxQperp(:,3),'k-s','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,FxQperp(:,4),'k-*','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
grid on;
legend('r=4','r=8','r=16','r=32','location','southwest');
xlabel('t');ylabel('$$\|\mathcal{P}_Q^\perp F^\star_X(t,Y,z)Q + X_d-Y\|_{M,\Xi}$$','interpreter','latex')
set(gca,'fontsize',18)
filename='../fig_output/residualadjointOtP.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 

figure(4)
semilogy(tvec,FxYperp(:,1),'k-o','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
hold on
semilogy(tvec,FxYperp(:,2),'k-x','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,FxYperp(:,3),'k-s','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
semilogy(tvec,FxYperp(:,4),'k-*','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
grid on;
legend('r=4','r=8','r=16','r=32','location','southwest');
xlabel('t');ylabel('$$\|\mathcal{P}_Y^\perp F^\star_X(t,Y,z)Q + X_d-Y\|_{M,\Xi}$$','interpreter','latex')
set(gca,'fontsize',18)
filename='../fig_output/residualadjointPtO.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 