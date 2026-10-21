close all;clear all;
load('ConvergenceGradients_Nh_961_N_500.mat');

figure(1)
loglog(dtvec,errOtP(1,:),'k--o','Linewidth',2,'Markersize',10);
hold on
loglog(dtvec,errOtP(2,:),'k--x','Linewidth',2,'Markersize',10);
loglog(dtvec,errOtP(3,:),'k--s','Linewidth',2,'Markersize',10);
loglog(dtvec,errOtP(4,:),'k--*','Linewidth',2,'Markersize',10);
loglog(dtvec,dtvec,'k','Linewidth',2);
grid on;xlabel('\Delta t');
ylabel('$$\|\nabla \widehat{J}(z)-g(z)\|_{L^2(I;L^2(\mathcal{D}))}$$','interpreter','latex'); 
set(gca,'fontsize',16);
legend('r=4','r=8','r=16','r=32','\Delta t','location','southeast');
filename='../fig_output/GradientOtP.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 

figure(2)
loglog(dtvec,errPtO(1,:),'k--o','Linewidth',2,'Markersize',10);
hold on
loglog(dtvec,errPtO(2,:),'k--x','Linewidth',2,'Markersize',10);
loglog(dtvec,errPtO(3,:),'k--s','Linewidth',2,'Markersize',10);
loglog(dtvec,errPtO(4,:),'k--*','Linewidth',2,'Markersize',10);
loglog(dtvec,dtvec,'k','Linewidth',2);
grid on;xlabel('\Delta t');
ylabel('$$\|\nabla \widehat{J}(z)-\nabla \widetilde{J}(z)\|_{L^2(I;L^2(\mathcal{D}))}$$','interpreter','latex'); 
set(gca,'fontsize',16);
legend('r=4','r=8','r=16','r=32','\Delta t','location','southeast');
filename='../fig_output/GradientPtO.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 

figure(3)
loglog(dtvec,errOtPAdjoint(1,:),'k--o','Linewidth',2,'Markersize',10);
hold on
loglog(dtvec,errOtPAdjoint(2,:),'k--x','Linewidth',2,'Markersize',10);
loglog(dtvec,errOtPAdjoint(3,:),'k--s','Linewidth',2,'Markersize',10);
loglog(dtvec,errOtPAdjoint(4,:),'k--*','Linewidth',2,'Markersize',10);
loglog(dtvec,dtvec,'k','Linewidth',2);
grid on;xlabel('\Delta t');
ylabel('$$\|P-Q^{OtP}\|_{L^2(I;L^2(\mathcal{D}))}$$','interpreter','latex'); 
set(gca,'fontsize',16);
legend('r=4','r=8','r=16','r=32','\Delta t','location','southeast');
filename='../fig_output/AdjointOtP.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 

figure(4)
loglog(dtvec,errPtOAdjoint(1,:),'k--o','Linewidth',2,'Markersize',10);
hold on
loglog(dtvec,errPtOAdjoint(2,:),'k--x','Linewidth',2,'Markersize',10);
loglog(dtvec,errPtOAdjoint(3,:),'k--s','Linewidth',2,'Markersize',10);
loglog(dtvec,errPtOAdjoint(4,:),'k--*','Linewidth',2,'Markersize',10);
loglog(dtvec,dtvec,'k','Linewidth',2);
grid on;xlabel('\Delta t');
ylabel('$$\|P-Q^{PtO}\|_{L^2(I;L^2(\mathcal{D}))}$$','interpreter','latex'); 
set(gca,'fontsize',16);
legend('r=4','r=8','r=16','r=32','\Delta t','location','southeast');
filename='../fig_output/AdjointPtO.eps';
exportgraphics(gcf,filename, 'ContentType', 'vector') 



    
    