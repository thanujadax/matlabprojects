% Figure - 14.16 and 14.17
% =========================================
% This script runs the global thresholding denoising experiment, and 
% demonstrates the true MSE as a function of the threshold T, and how 
% the SURE expression predicts well this value. The 
% measurement is created by taking the 'barbara' image and adding
% w.g.n. with sigma^2=400. This gives y. We aim to denoise by a 
% thresholding algorithm using an unecimated Haar with 3 levels of 
% resolution, 7:1 redundancy

% clear all; close all;

function []=Chapter_14_Global_SURE()

%================================================
%  S h o w i n g   t h e   s m o o t h e d   H T   c u r v e 
%================================================
x=-2:0.001:2; 
figure(1); clf; T=1; 
k=10; h=plot(x,x.*(x.^k)./(x.^k+T^k),'k:'); 
set(h,'LineWidth',2);
hold on; 
k=20; h=plot(x,x.*(x.^k)./(x.^k+T^k),'k--'); 
set(h,'LineWidth',2);
h=plot(x,((abs(x)>1).*x),'k'); 
set(h,'LineWidth',2);
grid on; 
legend({'k=10','k=20','Hard-Thresholding'},2); 
set(gca,'FontSize',14); 
% print -depsc2 Chapter_14_SmoothedHT.eps

%================================================
%  O r g a n i z i n g   t h e   D a t a   f o r   t h e   D e n o i s i n g
%================================================

y0=imread('barbara.png'); 
y0=double(y0(:));
n=length(y0); 

Haar=Generate_Haar_Matrix(sqrt(n));
[n,m]=size(Haar);
W=[0.25*ones(n*4,1); 0.5*ones(n*3,1)];

sigma=20; % noise power
y=y0+sigma*randn(n,1); % add noise
PSNRinput=10*log10(255^2/mean((y-y0).^2)); 

%================================================
%  D e n o i s i n g   a l g o r i t h m   a n d   r e s u l t s                      
%================================================

% First - The true and the SURE MSE
W2=W.^2; 
CurveT=zeros(1,201);
CurveSURE=zeros(1,201); 

 u=Haar'*y;
 u=u./W;
for T=0:1:200, 
    disp(T); pause(0.0001); 
    temp=(u/T).^20; 
    r=u.*(temp./(1+temp)); 
    yest=W.*r; 
    yest=Haar*yest; 
    CurveT(T+1)=sum((yest-y0).^2);  
    r=(temp.^2+21*temp)./(temp+1).^2; 
    CurveSURE(T+1)=yest'*yest-2*yest'*y+2*sigma^2*W2'*r; 
end;
CurveSURE=CurveSURE-CurveSURE(50)+CurveT(50);

figure(2); clf; 
h=plot(0:1:200,10*log10(255^2*length(y(:))./CurveSURE),'k');
set(h,'LineWidth',2);
hold on; 
h=plot(0:5:200,10*log10(255^2*length(y(:))./CurveT(1:5:end)),'ok');
set(h,'MarkerFaceColor','k');
h=xlabel('Threshold Value T');
set(h,'FontSize',14);
h=ylabel('PNSR');
set(h,'FontSize',14);
h=plot(0:1:200,ones(1,201)*PSNRinput,':k'); 
set(h,'LineWidth',2);
set(gca,'FontSize',14);
legend({'SURE PSNR','True PNSR','PSNR of y'}); 
axis([0 200 20 28]); 
% print -depsc2 Chapter_14_Global_T_SURE.eps

Topt=find(CurveSURE==min(CurveSURE),1)-1; 
u=Haar'*y;
u=u./W;
temp=(u/Topt).^20; 
r=u.*(temp./(1+temp)); 
yest=W.*r; 
yest=Haar*yest; 
PSNRoutput=10*log10(255^2/mean((yest-y0).^2)); 

disp(['PSNR input=',num2str(PSNRinput)]); 
disp(['PSNR output=',num2str(PSNRoutput)]); 

return;

%================================================
%================================================

function [Haar]=Generate_Haar_Matrix(n)

D1=sparse(n,n);
v=sparse([1 zeros(1,n-2), -1]);
for k=1:1:n
    D1(k,:)=v;
    v=[v(end),v(1:end-1)];
end;
D2=sparse(n,n);
v=[1 1 zeros(1,n-4), -1 -1];
for k=1:1:n
    D2(k,:)=v;
    v=[v(end),v(1:end-1)];
end;
S1=abs(D1);
S2=abs(D2);
Haar=[kron(S2,S2),kron(S2,D2),kron(D2,S2),kron(D2,D2),...
                             kron(S1,D1),kron(D1,S1),kron(D1,D1)];
for k=1:1:7*n^2
    Haar(:,k)=Haar(:,k)/sum(abs(Haar(:,k)));
end;
return;

