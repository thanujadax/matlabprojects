% Figures - 5.14, 5.15, and 5.16
% =========================================
% This program presents treates the inpainting problem locally as 
% the previous algorithm (Chapter_15_CoreInpainting2KSVD) , 
% but it operates on the image fingerprint.


function []=Chapter_15_CoreInpainting3KSVD()


% Gathering the data to train on from an image
% ========================================================

bb=8; % block size
K=256; % number of atoms in the dictionary
N=200; 

[IMin0,pp]=imread('fingerprint.png');
IMin0=im2double(IMin0(101:300,101:300));
IMin0=im2double(IMin0); 
IMin0 = IMin0*255; 

figure(1); clf; 
imagesc(IMin0); axis image; axis off; colormap(gray(256));
% print -depsc2 Chapter_15_CoreInpainting1.eps

% Adding noise
sigma=20;
IMin=IMin0+sigma*randn(N);

% Extracting the noisy patches
blkMatrixIm=im2col(IMin,[bb,bb],'sliding');

% ========================================================
disp('Inpainting results for the text mask:');
% ========================================================

disp('================================');

% The initial dictionary
Dict=zeros(bb,sqrt(K));
for k=0:1:sqrt(K)-1,
    V=cos([0:1:bb-1]'*k*pi/sqrt(K));
    if k>0, V=V-mean(V); end;
    Dict(:,k+1)=V/norm(V);
end;
Dict=kron(Dict,Dict);

% Creating the mask and extracting its patches
Mask=imread('TextMask.tif');
Mask=double(Mask(1:N,1:N)/255);
blkMask=im2col(Mask,[bb,bb],'sliding');

% Inpainting the Patches
Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 
Result1=RecoverImage(IMin,Dict,Coeff); 
Result1=max(min(Result1,255),0);
disp(['Recovery error =',num2str(sqrt(mean(mean((Result1-IMin0).^2))))]); 
figure(1); clf; 
imagesc([IMin0,IMin.*Mask; Result1, Result1]); 
axis image; colormap(gray(256)); drawnow;

for KSVDiter=1:1:15,    
    for atom=1:1:size(Dict,2)
        Omega=find(abs(Coeff(atom,:))>0);
        if isempty(Omega), continue; end; % this atom will be useless
        CoeffM=Coeff; 
        CoeffM(atom,:)=0; 
        Err=blkMask(:,Omega).*(blkMatrixIm(:,Omega)-Dict*CoeffM(:,Omega)); 
        for SVDiter=1:1:3
            dd=diag(1./(blkMask(:,Omega)*(Coeff(atom,Omega).^2)'))*...
                (Err*Coeff(atom,Omega)'); 
            Dict(:,atom)=dd/norm(dd); 
            Coeff(atom,Omega)=(dd'*Err)./sum((blkMask(:,Omega).*...
                     (dd*ones(1,length(Omega)))).^2,1); 
            % disp(norm(Err-blkMask(:,Omega).*(Dict(:,atom)*Coeff(atom,Omega)),'fro')); 
        end;
    end;
    Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,sigma*1.1);
    Result2=RecoverImage(IMin,Dict,Coeff);
    Result2=max(min(Result2,255),0);
    disp(['Iteration number: ',num2str(KSVDiter),...
            '     Recovery error =',num2str(sqrt(mean(mean((Result2-IMin0).^2))))]);    

end;

figure(1); clf; imagesc([IMin0, ones(200,5)*255, IMin.*Mask; ...
    ones(5,405)*255; Result1, ones(200,5)*255, Result2]);
axis image; colormap(gray(256)); axis off; drawnow;
% print -depsc2 Chapter_15_LocalInpaintingKSVDfngr1.eps

figure(2); image([abs(Result1-IMin0), ones(200,5)*1000, abs(IMin0-Result2)]*3);
axis image; colormap(gray(256)); axis off; drawnow;
% print -depsc2 Chapter_15_LocalInpaintingKSVDfngr2.eps

[b1,a1]=hist(IMin0(:)-Result1(:),-150:1:150);
[b2,a2]=hist(IMin0(:)-Result2(:),-150:1:150);
figure(3); 
h=plot(a1,b1); hold on; 
set(h,'LineWidth',2);
h=plot(a2,b2,'r--');
set(h,'LineWidth',2);
h=xlabel('Error');
set(h,'FontSize',14);
h=ylabel('Count');
set(h,'FontSize',14);
legend({'DCT error','K-SVD error'});
set(gca,'FontSize',14);
% print -depsc2 Chapter_15_LocalInpaintingKSVDfngr3.eps

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
