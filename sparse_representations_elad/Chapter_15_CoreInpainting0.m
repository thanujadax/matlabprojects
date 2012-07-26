% Figures - 5.10, 5.11, and 5.12
% =========================================
% This program presents the core idfea of inpainting using patch 
% processing. The treatment given here assumes no overlap between 
% the handled patches


function []=Chapter_15_CoreInpainting0()

% Gathering the data to train on from an image
% ========================================================

bb=8; % block size
K=256; % number of atoms in the dictionary

[IMin0,pp]=imread('peppers256.png');
IMin0=im2double(IMin0);
IMin0 = IMin0*255; 

figure(1); clf; 
imagesc(IMin0); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_CoreInpainting1.eps

% The dictionary
DCT=zeros(bb,sqrt(K));
for k=0:1:sqrt(K)-1,
    V=cos([0:1:bb-1]'*k*pi/sqrt(K));
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

% Adding noise
sigma=20;
IMin=IMin0+sigma*randn(256);

% Extracting the noisy patches
blkMatrixIm=im2col(IMin,[bb,bb],'distinct');
blkMatrixIm0=im2col(IMin0,[bb,bb],'distinct');

% ========================================================
% Part 1 - Generating a graph of the recovery quality versus %missing pixels
% ========================================================

Pos=randperm(256^2);
Error=zeros(100,1);
h=waitbar(0,'Testing with various missing percentages ...'); 
for percentage=5:5:100; 
    waitbar(percentage/100);
    % Generating the mask
    Mask=ones(256,256); 
    Mask(Pos(1:round(256*256*(100-percentage)/100)))=0; 
    blkMask=im2col(Mask,[bb,bb],'distinct');
    % Inpainting the Patches
    Coeff=OMPerrInpaint(DCT,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 
    % Computing the error 
    Reconst=min(max(DCT*Coeff,0),255); 
    Error(percentage)=sqrt(mean(mean((blkMatrixIm0-Reconst).^2))); 
    disp([percentage,Error(percentage)]); 
end;
close(h); 

figure(2); clf; 
h=plot(5:5:100,Error(5:5:100),'k');
set(h,'LineWidth',2);
hold on; 
h=plot(0:1:100,ones(1,101)*20,'--k');
set(h,'LineWidth',2);
axis([0 100 0 80]); 
h=xlabel('Percentage of remaining pixels');
set(h,'FontSize',14); 
h=ylabel('Root-Mean-Squared-Error'); 
set(h,'FontSize',14); 
set(gca,'FontSize',14); 
grid on;
legend({'Inpainting Results','Noise Level'}); 
% print -depsc2 Chapter_15_CoreInpaintingGRAPH.eps

% ========================================================
% Part 1 - Removing 25% of the pixels
% ========================================================

% Creating the mask and extracting its patches
Pos=randperm(256^2);
Mask=ones(256,256); 
Mask(Pos(1:256*64))=0; 
blkMask=im2col(Mask,[bb,bb],'distinct');

% Inpainting the Patches
Coeff=OMPerrInpaint(DCT,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 

% Creating the output image
Result1=col2im(DCT*Coeff,[bb,bb],[256,256],'distinct'); 
Result1=max(min(Result1,255),0);

figure(3); imagesc(IMin.*Mask); axis image; axis off; colormap(gray(256));
print -depsc2 Chapter_15_CoreInpainting2a.eps
figure(4); imagesc(Result1); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_CoreInpainting2b.eps

disp(['Recovery error =',num2str(sqrt(mean(mean((Result1-IMin0).^2))))]); 

% ========================================================
% Part 2 - Removing 50% of the pixels
% ========================================================

% Creating the mask and extracting its patches
Pos=randperm(256^2);
Mask=ones(256,256); 
Mask(Pos(1:256*128))=0; 
blkMask=im2col(Mask,[bb,bb],'distinct');

% Inpainting the Patches
Coeff=OMPerrInpaint(DCT,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 

% Creating the output image
Result2=col2im(DCT*Coeff,[bb,bb],[256,256],'distinct'); 
Result2=max(min(Result2,255),0);

figure(6); imagesc(IMin.*Mask); axis image; axis off; colormap(gray(256));
print -depsc2 Chapter_15_CoreInpainting3a.eps
figure(7); imagesc(Result2); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_CoreInpainting3b.eps

disp(['Recovery error =',num2str(sqrt(mean(mean((Result2-IMin0).^2))))]); 

% ========================================================
% Part 3 - Removing 75% of the pixels
% ========================================================

% Creating the mask and extracting its patches
Pos=randperm(256^2);
Mask=ones(256,256); 
Mask(Pos(1:256*192))=0; 
blkMask=im2col(Mask,[bb,bb],'distinct');

% Inpainting the Patches
Coeff=OMPerrInpaint(DCT,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 

% Creating the output image
Result3=col2im(DCT*Coeff,[bb,bb],[256,256],'distinct'); 
Result3=max(min(Result3,255),0);

figure(9); imagesc(IMin.*Mask); axis image; axis off; colormap(gray(256));
print -depsc2 Chapter_15_CoreInpainting4a.eps
figure(10); imagesc(Result3); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_CoreInpainting4b.eps

disp(['Recovery error =',num2str(sqrt(mean(mean((Result3-IMin0).^2))))]); 

return;

% ========================================================
% ========================================================

function [A]=OMPerrInpaint(D,X,Mask,errorGoal)
% ========================================================
% Sparse coding of a group of signals based on a given dictionary and specified representation
% error to get.
% input arguments: D - the dictionary
%                           X - the signals to represent
%                           errorGoal - the maximal allowed representation error
% output arguments: A - sparse coefficient matrix.
% ========================================================
[n,P]=size(X);
[n,K]=size(D);
E2 = errorGoal^2*n;
maxNumCoef = n/2;
A = sparse(size(D,2),size(X,2));
h=waitbar(0,'OMP on each example ...');
for k=1:1:P,
    waitbar(k/P);
    a=[];
    Mpos=find(Mask(:,k)); 
    E2M=E2*length(Mpos)/n;
    Dict=D(Mpos,:);
    W=1./sqrt(diag(Dict'*Dict)); 
    Dict=Dict*diag(W);
    x=X(Mpos,k);
    residual=x;
    indx = [];
    a = [];
    currResNorm2 = sum(residual.^2);
    j = 0;
    while currResNorm2>E2M && j < maxNumCoef,
        j = j+1;
        proj=Dict'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(Dict(:,indx(1:j)))*x;
        residual=x-Dict(:,indx(1:j))*a;
        currResNorm2=sum(residual.^2);
    end;
    if (~isempty(indx))
        A(indx,k)=a;
        A(:,k)=W.*A(:,k);
    end
end;
close(h); 
return;
