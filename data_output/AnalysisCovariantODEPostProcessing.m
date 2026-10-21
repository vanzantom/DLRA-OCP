close all;clear all;
filename=sprintf('AnalysisCovariantODE.mat');
load(filename);

loglog(dtvec,errfactorproj(1,:),'k-o','Linewidth',2,'Markersize',10);hold on;loglog(dtvec,errfactorfact(1,:),'k--o','Linewidth',2,'Markersize',10);
loglog(dtvec,errfactorproj(2,:),'k-x','Linewidth',2,'Markersize',10);loglog(dtvec,errfactorfact(2,:),'k--x','Linewidth',2,'Markersize',10);
loglog(dtvec,errfactorproj(3,:),'k-s','Linewidth',2,'Markersize',10);loglog(dtvec,errfactorfact(3,:),'k--s','Linewidth',2,'Markersize',10);
loglog(dtvec,errfactorproj(4,:),'k-*','Linewidth',2,'Markersize',10);loglog(dtvec,errfactorfact(4,:),'k--*','Linewidth',2,'Markersize',10);

grid on;
xlabel('\Delta t');
if local==true      
    loglog(dtvec,dtvec.^2,'k'); loglog(dtvec,dtvec.^5,'k--');
    legend(sprintf('rank %d',rvec(1)),'',sprintf('rank %d',rvec(2)),'',sprintf('rank %d',rvec(3)),'',sprintf('rank %d',rvec(4)),'','\Delta t^2','\Delta t^5','location','southeast');
    ylabel('Local truncation error');
    filename="../fig_output/TimeIntegrationCovariantODELTE.eps";

else
    loglog(dtvec,dtvec,'k'); loglog(dtvec,dtvec.^4,'k--');
    legend(sprintf('rank %d',rvec(1)),'',sprintf('rank %d',rvec(2)),'',sprintf('rank %d',rvec(3)),'',sprintf('rank %d',rvec(4)),'','\Delta t','\Delta t^4','location','southeast');
    ylabel('Error at final time'); 
    filename="../fig_output/TimeIntegrationCovariantODE.eps";

end
set(gca,'fontsize',16)
exportgraphics(gcf,filename, 'ContentType', 'vector') 

figure(2)
tvec=linspace(0,1,1/dtvec(end));
step=2^6;
semilogy(tvec,effectivecurvature(1,:),'k-o','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);hold on;
semilogy(tvec,effectivecurvature(2,:),'k-x','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10)
semilogy(tvec,effectivecurvature(3,:),'k-s','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10)
semilogy(tvec,effectivecurvature(4,:),'k-*','Linewidth',2,'MarkerIndices', 1:step:length(tvec),'Markersize',10);
grid on
ylim([0 50])
yticks([10 20 30 40 50])
legend(sprintf('rank %d',rvec(1)),sprintf('rank %d',rvec(2)),sprintf('rank %d',rvec(3)),sprintf('rank %d',rvec(4)));
xlabel('time');ylabel('Effective curvature');
set(gca,'fontsize',16)
filename="../fig_output/TimeIntegrationCovariantODEEffectiveCurvature.eps";
exportgraphics(gcf,filename, 'ContentType', 'vector') 