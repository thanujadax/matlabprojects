% Figures - 15.24, 15.25, and 15.26
% =========================================
% This program performs a simple salt-and-pepper denoising using 
% the inpainting idea


function []=Chapter_15_ImpulseNoiseRemoval()

% Gathering the data to train on from an image
% ========================================================

bb=8; % block size
K=256; % number of atoms in the dictionary
N=256; 

[IMin0,pp]=imread('peppers256.png');
% IMin0=im2double(IMin0(101:300,301:500));
IMin0=im2double(IMin0); 
IMin0 = IMin0*255; 

figure(1); clf; 
imagesc(IMin0); axis image; axis off; colormap(gray(256));

% Adding noise
portion=5; % percent noise of each type (salt & pepper)
pos=randperm(length(IMin0(:)));
Salt=pos(1:round(portion*length(IMin0(:))/100));
Pepper=pos(round(portion*length(IMin0(:))/100)+1:...
                   round(2*portion*length(IMin0(:))/100));
IMin=IMin0; 
IMin(Salt)=min(255,IMin(Salt)+50);
IMin(Pepper)=max(0,IMin(Pepper)-50); 

figure(1); clf; imagesc([IMin0, ones(N,5)*255, IMin]);
axis image; colormap(gray(256)); axis off; drawnow;
% print -depsc2 Chapter_15_ImpulseNoiseDCT0.eps

% Extracting the noisy patches
blkMatrixIm=im2col(IMin,[bb,bb],'sliding');

% The initial dictionary
Dict=zeros(bb,sqrt(K));
for k=0:1:sqrt(K)-1,
    V=cos([0:1:bb-1]'*k*pi/sqrt(K));
    if k>0, V=V-mean(V); end;
    Dict(:,k+1)=V/norm(V);
end;
Dict=kron(Dict,Dict);

% Estimating the mask and extracting its patches
EstimatedMask=abs(IMin-ordfilt2(IMin,5,ones(3)))<27; 
blkMask=im2col(EstimatedMask,[bb,bb],'sliding');

% Computing median filter all over the image
Result1=ordfilt2(IMin,5,ones(3)); 
disp(['Recovery error =',num2str(sqrt(mean(mean((Result1-IMin0).^2))))]); 

% Applying median only on the (estimated) masked pixels
Result2=Result1; 
Result2(EstimatedMask)=IMin(EstimatedMask);
disp(['Recovery error =',num2str(sqrt(mean(mean((Result2-IMin0).^2))))]); 

% Inpainting the Patches
Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,5); 
ResultTemp=RecoverImage(IMin,Dict,Coeff); 
ResultTemp=max(min(ResultTemp,255),0);
Result3=ResultTemp; 
Result3(EstimatedMask)=IMin(EstimatedMask);
disp(['Recovery error =',num2str(sqrt(mean(mean((Result3-IMin0).^2))))]); 

figure(1); clf; imagesc([IMin0, ones(N,5)*255, IMin; ...
                                  ones(5,2*N+5)*255; IMin0.*EstimatedMask, ones(N,5)*255, Result1; ...
                                  ones(5,2*N+5)*255; Result2, ones(N,5)*255, Result3]);
axis image; colormap(gray(256)); axis off; drawnow;
% print -depsc2 Chapter_15_ImpulseNoiseDCT1.eps

figure(2); image([abs(IMin-IMin0), ones(N,5)*1000, abs(IMin0-Result1); ...
                         ones(5,2*N+5)*1000; ...
                         abs(Result2-IMin0), ones(N,5)*1000, abs(IMin0-Result3)]*3);
axis image; colormap(gray(256)); axis off; drawnow;
% print -depsc2 Chapter_15_ImpulseNoiseDCT2.eps

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
