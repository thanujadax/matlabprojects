% Figures - 10.5 and 10.6
% =========================================
% This program plotes the rho function and the shrinkage that it induces

x=-1:0.0001:1; 
figure(1); clf; 
s=1; h=plot(x,abs(x)-s*log(1+abs(x)/s),'r');
set(h,'LineWidth',2);
hold on; 
s=0.1; h=plot(x,abs(x)-s*log(1+abs(x)/s),'g');
set(h,'LineWidth',2);
s=0.01; h=plot(x,abs(x)-s*log(1+abs(x)/s),'b');
set(h,'LineWidth',2);
s=0.001; h=plot(x,abs(x)-s*log(1+abs(x)/s),'m');
set(h,'LineWidth',2);
h=xlabel('x'); set(h,'FontSize',14); 
h=ylabel('\rho(x)');  set(h,'FontSize',14); 
set(gca,'Fontsize',14); 
legend({'s=1','s=0.1','s=0.01','s=0.001'}); 
% print -depsc2 Chapter_10_Shrinkage_Rho_Zibulevsky1.eps

figure(2); clf; 
s=1; h=plot(x,abs(x)-s*log(1+abs(x)/s),'r');
set(h,'LineWidth',2);
hold on; 
s=0.1; h=plot(x,abs(x)-s*log(1+abs(x)/s),'g');
set(h,'LineWidth',2);
s=0.01; h=plot(x,abs(x)-s*log(1+abs(x)/s),'b');
set(h,'LineWidth',2);
s=0.001; h=plot(x,abs(x)-s*log(1+abs(x)/s),'m');
set(h,'LineWidth',2);
axis([-0.05 0.05 0 0.05]); 
h=xlabel('x'); set(h,'FontSize',14); 
h=ylabel('\rho(x)');  set(h,'FontSize',14); 
set(gca,'Fontsize',14); 
legend({'s=1','s=0.1','s=0.01','s=0.001'}); 
% print -depsc2 Chapter_10_Shrinkage_Rho_Zibulevsky2.eps

figure(3); clf; 
x=-3:0.0001:3; 
lambda=1; 
s=1; 
Res = (x>=0).*(x-s-lambda+sqrt((x-s-lambda).^2+4*x*s))/2 ...
          -(x<0).*(-x-s-lambda+sqrt((-x-s-lambda).^2-4*x*s))/2; 
h=plot(x,Res,'r'); 
set(h,'LineWidth',2);
hold on;
s=0.1; 
Res = (x>=0).*(x-s-lambda+sqrt((x-s-lambda).^2+4*x*s))/2 ...
          -(x<0).*(-x-s-lambda+sqrt((-x-s-lambda).^2-4*x*s))/2; 
h=plot(x,Res,'g'); 
set(h,'LineWidth',2);
s=0.01; 
Res = (x>=0).*(x-s-lambda+sqrt((x-s-lambda).^2+4*x*s))/2 ...
          -(x<0).*(-x-s-lambda+sqrt((-x-s-lambda).^2-4*x*s))/2; 
h=plot(x,Res,'b'); 
set(h,'LineWidth',2);
s=0.001; 
Res = (x>=0).*(x-s-lambda+sqrt((x-s-lambda).^2+4*x*s))/2 ...
          -(x<0).*(-x-s-lambda+sqrt((-x-s-lambda).^2-4*x*s))/2; 
h=plot(x,Res,'m'); 
set(h,'LineWidth',2);
h=plot(-3:0.1:3,-3:0.1:3,':');
set(h,'LineWidth',2);
h=xlabel('x'); set(h,'FontSize',14); 
h=ylabel('S_{\rho,\lambda}(x)');  set(h,'FontSize',14); 
set(gca,'Fontsize',14); 
legend({'s=1','s=0.1','s=0.01','s=0.001','y=x'},2); 
% print -depsc2 Chapter_10_Shrinkage_Rho_Zibulevsky3.eps

