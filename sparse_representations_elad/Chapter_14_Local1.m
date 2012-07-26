% Figure - 14.3, 14.4, and 14.5
% =========================================
% This script runs an image denoising experiment with patches. The 
% algorithm considered applies the unitary 2D-DCt on each patch, and 
% trains a shrinkage curve per each coefficient. The training applied 
% here is fitting patches (input to output (disregarding the overlaps)


clear all; close all;

n=6; % patch size

%================================================
%  O r g a n i z i n g   t h e   T r a i n i n g   D a t a
%================================================

y0=imread('lena.png'); 
y0=double(y0(101:300,101:300));
% y0=imread('barbara.png'); 
% y0=double(y0);
N=size(y0,1);

figure(1); clf; 
image(reshape(y0,[N,N]));
colormap(gray(256)); axis image; axis off; 
% print -depsc2 Chapter_14_Local1_Train.eps

noise=randn(N,N);
sigma=20; % noise power
y=y0+sigma*noise; % add noise

%================================================
% L e a r n   t h e   S h r i n k a g e   c u r v e s   b y   p a t c h e s   d i r e c t l y 
%                                T h e   U N I T A R Y   c a s e
%================================================

% accumulate the clean & noisy n*n patches
P0=zeros(n^2,(N-n+1)^2);
P=zeros(n^2,(N-n+1)^2);
count=1; 
for j=1:1:N-n+1
    for i=1:1:N-n+1
        patch=y0(i:i+n-1,j:j+n-1);
        P0(:,count)=patch(:);
        patch=y(i:i+n-1,j:j+n-1);
        P(:,count)=patch(:); 
        count=count+1;
    end;
end;
P=(P-127)/128;
P0=(P0-127)/128;

% the unitary transform to apply
T=kron(dctmtx(n),dctmtx(n));
P=T*P; 
P0=T*P0; 

% Fitting the LUT curves
figure(1); clf; 
Coef=zeros(4,n^2); 
Limit=zeros(n^2,2);
for k=1:1:n^2,
    Limit(k,:)=[min(P(k,:)),max(P(k,:))];
    mp=max(abs(P(k,:)));
    Vander=[P(k,:).^1; P(k,:).^3; P(k,:).^5; P(k,:).^7]'; 
    Coef(:,k)=inv(Vander'*Vander)*(Vander'*P0(k,:)'); 
    v=(-mp:0.01:mp)'; 
    subplot(n,n,k); 
    plot(v,[v.^1, v.^3, v.^5, v.^7]*Coef(:,k),'r');
    hold on; plot([-mp mp],[-mp mp],':');
    axis equal; 
    axis([-mp ,mp -mp mp]); 
    drawnow;
end;
% print -depsc2 Chapter_14_Local1_Curves.eps

%================================================
%  F i l t e r   a   n e w   i m a g e   ( B a r b a r a )
%================================================

y0=imread('barbara.png'); 
y0=double(y0);
N=size(y0,1);

noise=randn(N,N);
sigma=20; % noise power
y=y0+sigma*noise; % add noise
PSNRinput=10*log10(255^2/mean((y(:)-y0(:)).^2)); 

% accumulate the noisy n*n patches
P=zeros(n^2,(N-n+1)^2);
count=1; 
for j=1:1:N-n+1
    for i=1:1:N-n+1
        patch=y(i:i+n-1,j:j+n-1);
        P(:,count)=patch(:); 
        count=count+1;
    end;
end;
P=(P-127)/128;
P=T*P; 

% Denoising the image
P0hat=P; 
for k=1:1:n^2,
    Mask=(P(k,:)>=Limit(k,1) & P(k,:)<=Limit(k,2)); 
    Vander=[P(k,:).^1; P(k,:).^3; P(k,:).^5; P(k,:).^7]'; 
    P0hat(k,:)=Vander*Coef(:,k);
    P0hat(k,:)=P0hat(k,:).*Mask+(1-Mask).*P(k,:);
end;
P0hat=T'*P0hat;
P0hat=128*P0hat+127; 

yout=zeros(N,N); 
Weight=zeros(N,N); 
i=1; j=1;
for k=1:1:(N-n+1)^2,
    patch=reshape(P0hat(:,k),[n,n]); 
    yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch; 
    Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1; 
    if i<N-n+1 
        i=i+1; 
    else
        i=1; j=j+1; 
    end;
end;
yout=yout./Weight; 

PSNRoutput=10*log10(255^2/mean((yout(:)-y0(:)).^2)); 
disp([PSNRinput,PSNRoutput]);

figure(2); clf; 
image(reshape(y0,[512,512]));
colormap(gray(256)); axis image; axis off; 
% print -depsc2 Chapter_14_Local1_Orig.eps

figure(3); clf; 
image(reshape(y,[512,512]));
colormap(gray(256)); axis image; axis off; 
% print -depsc2 Chapter_14_Local1_Noisy.eps

figure(4); clf; 
image(reshape(yout,[512,512]));
colormap(gray(256)); axis image; axis off; 
% print -depsc2 Chapter_14_Local1_Result.eps


