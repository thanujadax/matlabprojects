% Figures - NONE
% =========================================
% This program presents treates the inpainting problem locally as 
% the two previous algorithms, but introduces a leearned dictionary. 
% The results built here are used to fill Table 15.2


function []=Chapter_15_CoreInpainting2KSVD()

% Gathering the data to train on from an image
% ========================================================

bb=8; % block size
K=256; % number of atoms in the dictionary
N=256; 

[IMin0,pp]=imread('peppers256.png');
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
% Part 1 - Removing 25% of the pixels
% ========================================================

disp('Inpainting results: 25 missing pixels');
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
Pos=randperm(N^2);
Mask=ones(N,N); 
Mask(Pos(1:N^2/4))=0; 
blkMask=im2col(Mask,[bb,bb],'sliding');

% Inpainting the Patches
Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 
Result1=RecoverImage(IMin,Dict,Coeff); 
Result1=max(min(Result1,255),0);
disp(['Recovery error =',num2str(sqrt(mean(mean((Result1-IMin0).^2))))]); 

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
            dd=dd/norm(dd); 
            Dict(:,atom)=dd; 
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

% ========================================================
% Part 2 - Removing 50% of the pixels
% ========================================================

disp('Inpainting results: 50 missing pixels');
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
Pos=randperm(N^2);
Mask=ones(N,N); 
Mask(Pos(1:N^2/2))=0; 
blkMask=im2col(Mask,[bb,bb],'sliding');

% Inpainting the Patches
Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 
Result3=RecoverImage(IMin,Dict,Coeff); 
Result3=max(min(Result3,255),0);
disp(['Recovery error =',num2str(sqrt(mean(mean((Result3-IMin0).^2))))]); 

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
            dd=dd/norm(dd); 
            Dict(:,atom)=dd; 
            Coeff(atom,Omega)=(dd'*Err)./sum((blkMask(:,Omega).*...
                     (dd*ones(1,length(Omega)))).^2,1); 
            % disp(norm(Err-blkMask(:,Omega).*(Dict(:,atom)*Coeff(atom,Omega)),'fro')); 
        end;
    end;
    Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,sigma*1.1);
    Result4=RecoverImage(IMin,Dict,Coeff);
    Result4=max(min(Result4,255),0);
    disp(['Iteration number: ',num2str(KSVDiter),...
            '     Recovery error =',num2str(sqrt(mean(mean((Result4-IMin0).^2))))]);
end;

% ========================================================
% Part 3 - Removing 75% of the pixels
% ========================================================

disp('Inpainting results: 75 missing pixels');
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
Pos=randperm(N^2);
Mask=ones(N,N); 
Mask(Pos(1:3*N^2/4))=0; 
blkMask=im2col(Mask,[bb,bb],'sliding');

% Inpainting the Patches
Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,sigma*1.1); 
Result5=RecoverImage(IMin,Dict,Coeff); 
Result5=max(min(Result5,255),0);
disp(['Recovery error =',num2str(sqrt(mean(mean((Result5-IMin0).^2))))]); 

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
            dd=dd/norm(dd); 
            Dict(:,atom)=dd; 
            Coeff(atom,Omega)=(dd'*Err)./sum((blkMask(:,Omega).*...
                     (dd*ones(1,length(Omega)))).^2,1); 
            % disp(norm(Err-blkMask(:,Omega).*(Dict(:,atom)*Coeff(atom,Omega)),'fro')); 
        end;
    end;
    Coeff=OMPerrInpaint(Dict,blkMatrixIm.*blkMask,blkMask,sigma*1.1);
    Result6=RecoverImage(IMin,Dict,Coeff);
    Result6=max(min(Result6,255),0);
    disp(['Iteration number: ',num2str(KSVDiter),...
            '     Recovery error =',num2str(sqrt(mean(mean((Result6-IMin0).^2))))]);
end;

save temp.mat; 

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
