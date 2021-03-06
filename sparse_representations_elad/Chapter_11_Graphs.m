% Figures - 11.1, 11.2, 11.4, 11.5, and then 11.6, 11.7, 11.8, 11.9
% =========================================
% This program tests the Oreacle, MAP, MMSE and their approximations
% for a very low-dimensional case, where the exhaustive expressions
% can be evaluated. Then it turns to test a higher-dimensional case
% where only the approximations are relevant. Finally, it shows the 
% Random-OMP behavior

% =========================================
% Test 1 - very low dimensional case, in which all the estimators 
%             can be applied 
% =========================================

Results=Chapter_11_Synthetic(20,30,3,1000,25,ones(1,10));
save Results1
V=(0.1:0.1:2)';

% Figure 1 - The Oracle
figure(1); clf;
h=plot(V,mean(Results(:,:,2),2)./(V.^2*20),'k');
set(h,'LineWidth',2);
hold on;
h=plot(V,mean(Results(:,:,3),2)./(V.^2*20),'k:');
set(h,'LineWidth',2);
h=plot(V,3/20*ones(20,1),'k--'); 
set(h,'LineWidth',2);
h=xlabel('\sigma_e');
set(h,'FontSize',14);
h=ylabel('Relative Mean-Squared-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
legend({'Oracle (empirical)','Oracle (formula)','k/n factor'});
axis([0.1 2 0 0.2]);
% print -depsc2 Chapter_11_Oracle.eps

% Figure 2 - The MAP and its approximation
figure(2); clf;
h=plot(V,mean(Results(:,:,2),2)./(V.^2*20),'k:');
set(h,'LineWidth',2);
hold on;
h=plot(V,mean(Results(:,:,4),2)./(V.^2*20),'k');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,6),2)./(V.^2*20),'ko');
set(h,'LineWidth',2);
set(h,'MarkerFaceColor','k');
h=plot(V,3/20*ones(20,1),'k--'); 
set(h,'LineWidth',2);
h=xlabel('\sigma_e');
set(h,'FontSize',14);
h=ylabel('Relative Mean-Squared-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
legend({'Oracle (empirical)','Exact MAP (empirical)',...
            'Approx. MAP by OMP','k/n factor'});
axis([0.1 2 0 0.6]);
% print -depsc2 Chapter_11_MAP.eps

% Figure 3 - The MMSE and its approximation
figure(3); clf;
h=plot(V,mean(Results(:,:,2),2)./(V.^2*20),'k:');
set(h,'LineWidth',2);
hold on;
h=plot(V,mean(Results(:,:,4),2)./(V.^2*20),'k--');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,7),2)./(V.^2*20),'k');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,9),2)./(V.^2*20),'ko');
set(h,'MarkerFaceColor','k');
plot(V,3/20*ones(20,1),'k-.'); 
h=xlabel('\sigma_e');
set(h,'FontSize',14);
h=ylabel('Relative Mean-Squared-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
legend({'Oracle (empirical)','Exact MAP (empirical)',...
            'Exact MMSE (empirical)','Approx. MMSE - RandOMP','k/n factor'});
axis([0.1 2 0 0.6]);
% print -depsc2 Chapter_11_MMSE.eps

% Figure 4 - The MAP and MMSE formulas
figure(4); clf;
h=plot(0.1:0.1:2,mean(Results(:,:,2),2)./(V.^2*20),'k:');
set(h,'LineWidth',2);
hold on;
h=plot(V,mean(Results(:,:,4),2)./(V.^2*20),'k--');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,5),2)./(V.^2*20),'k+');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,7),2)./(V.^2*20),'k');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,8),2)./(V.^2*20),'kx');
set(h,'LineWidth',2);
h=plot(V,3/20*ones(20,1),'k-.'); 
set(h,'LineWidth',2);
h=xlabel('\sigma_e');
set(h,'FontSize',14);
h=ylabel('Relative Mean-Squared-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
axis([0.1 2 0 0.6]);
legend({'Oracle (empirical)','Exact MAP (empirical)','Exact MAP (formula)',...
                   'Exact MMSE (empirical)','Exact MMSE (formula)','k/n factor'});
% print -depsc2 Chapter_11_Formula.eps

% ================================================
% Test 2 - A higher dimensional case, in which only approximations (and the 
%             oracle) can be tested 
% ================================================

% Figure 5 - Practical Estimators
Results=Chapter_11_Synthetic(200,400,20,100,25,[1,1,0,0,0,1,0,0,1,1]);

figure(5); clf;
h=plot(V,mean(Results(:,:,2),2)./(V.^2*200),'b');
set(h,'LineWidth',2);
hold on;
h=plot(V,mean(Results(:,:,6),2)./(V.^2*200),'g');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,9),2)./(V.^2*200),'r');
set(h,'LineWidth',2);
h=plot(V,mean(Results(:,:,10),2)./(V.^2*200),'c');
set(h,'LineWidth',2);
h=plot(V,20/200*ones(20,1),'k-.'); 
set(h,'LineWidth',2);
legend({'Oracle (empirical)','Approx. MAP by OMP',...
                   'Approx. MMSE by RandOMP',...
                   'Sparsified RandOMP','k/n factor'});
h=xlabel('\sigma_e');
set(h,'FontSize',14);
h=ylabel('Relative Mean-Squared-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
axis([0.1 2 0 0.5]);
% print -depsc2 Chapter_11_Practical.eps

% Figure 6 - Effect of Random-OMP runs
Results1=Chapter_11_Synthetic(200,400,20,100,1,[1,0,0,0,0,0,0,0,1,1]);
Results2=Chapter_11_Synthetic(200,400,20,100,2,[1,0,0,0,0,0,0,0,1,1]);
Results3=Chapter_11_Synthetic(200,400,20,100,4,[1,0,0,0,0,0,0,0,1,1]);
Results4=Chapter_11_Synthetic(200,400,20,100,8,[1,0,0,0,0,0,0,0,1,1]);
Results5=Chapter_11_Synthetic(200,400,20,100,16,[1,0,0,0,0,0,0,0,1,1]);
Results6=Chapter_11_Synthetic(200,400,20,100,32,[1,0,0,0,0,0,0,0,1,1]);
save Results2

figure(6); clf;
h=plot(V,mean(Results1(:,:,9),2)./(V.^2*200),'k--');
set(h,'LineWidth',2);
hold on;
h=plot(V,mean(Results2(:,:,9),2)./(V.^2*200),'b--');
set(h,'LineWidth',2);
h=plot(V,mean(Results3(:,:,9),2)./(V.^2*200),'c');
set(h,'LineWidth',2);
h=plot(V,mean(Results4(:,:,9),2)./(V.^2*200),'m');
set(h,'LineWidth',2);
h=plot(V,mean(Results5(:,:,9),2)./(V.^2*200),'r');
set(h,'LineWidth',2);
h=plot(V,mean(Results6(:,:,9),2)./(V.^2*200),'k');
set(h,'LineWidth',2);
legend({'J=1','J=2','J=4','J=8','J=16','J=32'});
h=xlabel('\sigma_e');
set(h,'FontSize',14);
h=ylabel('Random-OMP Relative Mean-Squared-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
axis([0.1 2 0 0.4]);
% print -depsc2 Chapter_11_RandOMP_runs.eps

% ================================================
% Test 3 - is random-OMP sparse?  
% ================================================

[A,S0,x0,y,x_randomp,x_sp]=Chapter_11_RandomOMP(200,400,10,100,0.2);
SNR=(A*x0)'*(A*x0)/((A*x0-y)'*(A*x0-y));
Error=sqrt((x0-x_randomp)'*(x0-x_randomp)/length(x0));
disp([SNR,Error/0.2]);

figure(7); clf; 
h=plot(x_randomp,'k');
set(h,'LineWidth',1);
hold on;
h=plot(S0,x0(S0),'ko');
set(h,'MarkerfaceColor','k');
legend({'Random-OMP Result','x_0'});
h=xlabel('index');
set(h,'FontSize',14);
h=ylabel('Amplitude');
set(h,'FontSize',14);
set(gca,'FontSize',14);
% print -depsc2 Chapter_11_RandOMP_outcome1.eps

figure(8); clf; 
h=semilogy(sort(abs(x0)+1e-10,'descend'),'k');
set(h,'LineWidth',2);
hold on;
h=semilogy(sort(abs(x_randomp)+1e-10,'descend'),'--k');
set(h,'LineWidth',2);
legend({'x_0','Random-OMP Result'});
h=xlabel('Sorted index');
set(h,'FontSize',14);
h=ylabel('Absolute Amplitude');
set(h,'FontSize',14);
set(gca,'FontSize',14);
axis([1 400 1e-10 2]);
grid on;
% print -depsc2 Chapter_11_RandOMP_outcome2.eps


