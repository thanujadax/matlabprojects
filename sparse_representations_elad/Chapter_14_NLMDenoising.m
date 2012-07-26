% Figure - 14.14
% =========================================
% This script offers a crude way to perform the NLM filter on Barbara.
% Requirement: be patient (or write your own code)



n=3; % patch size is (2n+1)^2
sigma=20; % noise power
M=8; % neighborhood size (2M+1)^2
T=sigma^2; % scale of the L2 distances

% Gather the data from the noisy Barbara
y0=imread('barbara.png'); 
y0=double(y0);
N=size(y0,1);
noise=randn(N,N);
y=y0+sigma*noise; % add noise
PSNRinput=10*log10(255^2/mean((y(:)-y0(:)).^2)); 

% Denoise by NLM
yout=zeros(N,N);
h=waitbar(0,'NLM filtering ...');
y=[y(:,M+n:-1:1), y, y(:,end-M-n+1:end)];
y=[y(M+n:-1:1,:); y; y(end-M-n+1:end,:)];
for i=n+M+1:1:N+M+n
    waitbar((i-M-n)/N);
    for j=n+M+1:1:N+M+n
        Center=y(i-n:i+n,j-n:j+n); 
        Weights=zeros(2*M+1,2*M+1); 
        for p=-M:M
            for q=-M:M
                Patch=y(i+p-n:i+p+n,j+q-n:j+q+n);
                dist2=mean((Patch(:)-Center(:)).^2); 
                Weights(p+M+1,q+M+1)=exp(-dist2/T); 
            end;
        end;
        Weights=Weights/sum(Weights(:)); 
        yout(i-n-M,j-n-M)=sum(sum(y(i-M:i+M,j-M:j+M).*Weights)); 
    end;
end;
close(h); 

y=y(M+n+1:M+n+N,M+n+1:M+n+N); 
PSNRoutput=10*log10(255^2/mean((yout(:)-y0(:)).^2)); 
disp([PSNRinput,PSNRoutput]);
figure(1); clf; 
image([y0,y,yout]); colormap(gray(256)); axis image; axis off; drawnow;
figure(2); clf; 
image(yout); colormap(gray(256)); axis image; axis off; drawnow;
% print -depsc2 Chapter_14_NLMres.eps


