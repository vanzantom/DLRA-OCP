close all;clear all;
load('ConvergenceOptimization_Nh_961_N_500.mat');

maxit=100;
step=10;
figure(1)
Jopt=Full.Jvec(end);
semilogy((0:maxit),(Full.Jvec(1:maxit+1)-Jopt)/Jopt,'k-d','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
hold on
semilogy((0:OtP.iter),(OtP.JvecRef(1,:)-Jopt)/Jopt,'k-o','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
semilogy((0:OtP.iter),(OtP.JvecRef(2,:)-Jopt)/Jopt,'k-x','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
semilogy((0:OtP.iter),(OtP.Jvec(3,:)-Jopt)/Jopt,'k-s','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
semilogy((0:OtP.iter),(OtP.Jvec(4,:)-Jopt)/Jopt,'k-*','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
grid on;xlabel('Iterations'); ylabel('$$|\widehat{J}(z^k)-\widehat{J}^\star|/\widehat{J}^\star$$','interpreter','latex');
legend('full rank','r=4','r=8','r=16','r=32','location','northeast');title('Optimize-then-Project')
set(gca,'fontsize',16)
filename='../fig_output/ConvergenceOptimizationOtP.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 

maxit=100;
figure(2)
Jopt=Full.Jvec(end);
semilogy((0:maxit),(Full.Jvec(1:maxit+1)-Jopt)/Jopt,'k-d','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
hold on
semilogy((0:OtP.iter),(PtO.JvecRef(1,:)-Jopt)/Jopt,'k-o','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
semilogy((0:OtP.iter),(PtO.JvecRef(2,:)-Jopt)/Jopt,'k-x','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
semilogy((0:OtP.iter),(PtO.Jvec(3,:)-Jopt)/Jopt,'k-s','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
semilogy((0:OtP.iter),(PtO.Jvec(4,:)-Jopt)/Jopt,'k-*','Linewidth',2,'MarkerIndices', 1:step:maxit,'Markersize',10);
grid on;xlabel('Iterations'); ylabel('$$|\widehat{J}(z^k)-\widehat{J}^\star|/\widehat{J}^\star$$','interpreter','latex');
legend('full rank','r=4','r=8','r=16','r=32','location','northeast');
title('Project-then-Optimize')
set(gca,'fontsize',16)
filename='../fig_output/ConvergenceOptimizationPtO.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 
