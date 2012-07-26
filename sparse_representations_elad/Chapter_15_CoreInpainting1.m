% Figures - 5.13
% =========================================
% This program presents the core idfea of inpainting using patch 
% processing. The treatment given here assumes the existance of
% overlaps between the handled patches, averaging their values.


function []=Chapter_15_CoreInpainting1()

% Gathering the data to train on from an image
% ========================================================

bb=8; % block size
K=256; % number of atoms in the dictionary

[IMin0,pp]=imread('peppers256.png');
IMin0=im2double(IMin0);
IMin0 = IMin0*255; 

% figure(1); clf; 
% imagesc(IMin0); axis image; axis off; colormap(gray(256));
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
blkMatrixIm=im2col(IMin,[bb,bb],'sliding');

% ========================================================
% Part 1 - Removing 25% of the pixels
% ========================================================

% Creating the mask and extracting its patches
Pos=randperm(256^2);
Mask=ones(256,256); 
Mask(Pos(1:256*64))=0; 
blkMask=im2col(Mask,[bb,bb],'sliding');

% Inpainting the Patches
Coeff=OMPerrInpaint(DCT,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 

% Creating the output image
Result1=RecoverImage(IMin,DCT,Coeff); 
Result1=max(min(Result1,255),0);
figure(2); imagesc([IMin0,IMin;IMin.*Mask,Result1]); 
axis image; axis off; colormap(gray(256))

figure(3); imagesc(IMin.*Mask); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_LocalInpaintingDCT1a.eps
figure(4); imagesc(Result1); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_LocalInpaintingDCT1b.eps

disp(['Recovery error =',num2str(sqrt(mean(mean((Result1-IMin0).^2))))]); 

% ========================================================
% Part 2 - Removing 50% of the pixels
% ========================================================

% Creating the mask and extracting its patches
Pos=randperm(256^2);
Mask=ones(256,256); 
Mask(Pos(1:256*128))=0; 
blkMask=im2col(Mask,[bb,bb],'sliding');

% Inpainting the Patches
Coeff=OMPerrInpaint(DCT,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 

% Creating the output image
Result2=RecoverImage(IMin,DCT,Coeff); 
Result2=max(min(Result2,255),0);
figure(5); imagesc([IMin0,IMin;IMin.*Mask,Result2]); 
axis image; axis off; colormap(gray(256))

figure(6); imagesc(IMin.*Mask); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_LocalInpaintingDCT2a.eps
figure(7); imagesc(Result2); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_LocalInpaintingDCT2b.eps

disp(['Recovery error =',num2str(sqrt(mean(mean((Result2-IMin0).^2))))]); 

% ========================================================
% Part 3 - Removing 75% of the pixels
% ========================================================

% Creating the mask and extracting its patches
Pos=randperm(256^2);
Mask=ones(256,256); 
Mask(Pos(1:256*192))=0; 
blkMask=im2col(Mask,[bb,bb],'sliding');

% Inpainting the Patches
Coeff=OMPerrInpaint(DCT,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 

% Creating the output image
Result3=RecoverImage(IMin,DCT,Coeff); 
Result3=max(min(Result3,255),0);
figure(8); imagesc([IMin0,IMin;IMin.*Mask,Result3]); 
axis image; axis off; colormap(gray(256))

figure(9); imagesc(IMin.*Mask); axis image; axis off; colormap(gray(256));
print -depsc2 Chapter_15_LocalInpaintingDCT3a.eps
figure(10); imagesc(Result3); axis image; axis off; colormap(gray(256));
print -depsc2 Chapter_15_LocalInpaintingDCT3b.eps

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

% ========================================================
% ========================================================

function [yout]=RecoverImage(y,D,CoefMatrix)
% ========================================================
% ========================================================
N=size(y,1); 
n=sqrt(size(D,1)); 
yout=zeros(N,N); 
Weight=zeros(N,N); 
i=1; j=1;
for k=1:1:(N-n+1)^2,
    patch=reshape(D*CoefMatrix(:,k),[n,n]); 
    yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch; 
    Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1; 
    if i<N-n+1 
        i=i+1; 
    else
        i=1; j=j+1; 
    end;
end;
yout=yout./Weight; 
return;
