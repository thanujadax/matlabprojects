% Figures - 14.1 and 14.2
% =========================================
% This script runs the denoising experiment. The details of the 
% experiment are the following:
% 1. Measurement created by taking the 'barbara' image and adding
%     w.g.n. with sigma^2=100. This gives y.
% 2. We aim to denoise by a thresholding algorithm using an 
%     unecimated Haar with 3 levels of resolution, 7:1 redundancy

% clear all; close all;

function []=Chapter_14_Global()

%================================================
%  O r g a n i z i n g   t h e   D a t a
%================================================

% --------------------------------- Part 1 - creation of data ------------------------------------

y0=imread('barbara.png'); 
% y0=y0(1:100,1:100);
y0=double(y0(:));
n=length(y0); 

Haar=Generate_Haar_Matrix(sqrt(n));
[n,m]=size(Haar);
W=[0.25*ones(n*4,1); 0.5*ones(n*3,1)];

noise=randn(n,1);
sigma=20; % noise power
y=y0+sigma*noise; % add noise
PSNRinput=10*log10(255^2/mean((y-y0).^2)); 

%================================================
%  D e n o i s i n g   a l g o r i t h m   a n d   r e s u l t s                      
%================================================

MM=spdiags(W,[0],m,m);
WW=spdiags(1./W,[0],m,m);
CurveT=zeros(1,200);
for T=0:1:200, 
    yest=Haar*((abs(Haar'*y)>T*W).*(Haar'*y));
    % yest=(abs(WW*Haar'*y)>T).*(WW*Haar'*y); 
    % yest=MM*yest; 
    % yest=Haar*yest;
    CurveT(T+1)=10*log10(255^2/mean((yest-y0).^2)); 
    disp(CurveT(T+1)); 
end;

figure(1); clf; 
h=plot(0:1:200,CurveT,'k');
set(h,'LineWidth',2);
h=xlabel('Threshold Value T');
set(h,'FontSize',14);
h=ylabel('PNSR');
set(h,'FontSize',14);
hold on;
h=plot(0:1:200,ones(1,201)*PSNRinput,'k:'); 
set(h,'LineWidth',2);
set(gca,'FontSize',14);
legend({'PNSR for various T values','PSNR of y'}); 
axis([0 200 20 28]); 
% print -depsc2 Chapter_14_Global_Tchoice.eps

Topt=find(CurveT==max(CurveT),1)-1; 
yest=Haar*((abs(Haar'*y)>Topt*W).*(Haar'*y)); 
PSNRoutput=10*log10(255^2/mean((yest-y0).^2)); 

figure(1); clf; 
image(reshape(y0,[512,512]));
colormap(gray(256)); axis image; axis off; 
% print -depsc2 Chapter_14_Global_Orig.eps

figure(2); clf; 
image(reshape(y,[512,512]));
colormap(gray(256)); axis image; axis off; 
% print -depsc2 Chapter_14_Global_Noisy.eps

figure(2); clf; 
image(reshape(yest,[512,512]));
colormap(gray(256)); axis image; axis off; 
% print -depsc2 Chapter_14_Global_Result.eps

disp(['PSNR input=',num2str(PSNRinput)]); 
disp(['PSNR output=',num2str(PSNRoutput)]); 
disp(['Sparsity of the representation: ',...
                   num2str(nnz((abs(Haar'*y)>Topt*W))/m)]);

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

