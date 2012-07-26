function []=Chapter_15_ScaleUp_Training_and_Testing2(ImageIN)

%============================================
%                              Training the two dictionaries
%============================================

% Setting parameters
n=9; % block size
m=1000; % number of atoms in the dictionary
s=2; % scale-down factor
dd=3; % margins in the image to avoid (dd*s to each side)
L=3; % number of atoms to use in the representation

% Preparing the low-and high resolution images
Yh=imread(ImageIN);
N=size(Yh,1);
N=floor(N/s^2)*s^2; 
Yh=im2double(Yh(1:N,1:N)); % so that it scales down to an integer size
Yh=Yh*255;

% Creating the low-resolution image to train on
Zl=conv2(Yh,[1 2 1]/4,'same');
Zl=conv2(Zl,[1 2 1]'/4,'same');
Yh=Zl(1:s:end,1:s:end); 

% Creating the low-low-resolution image 
Zl=conv2(Yh,[1 2 1]/4,'same');
Zl=conv2(Zl,[1 2 1]'/4,'same');
Zl=Zl(1:s:end,1:s:end); 

% Upscaling Zl to the original resolution
N=size(Yh,1);
[posY,posX]=meshgrid(1:s:N,1:s:N); 
[posY0,posX0]=meshgrid(1:1:N,1:1:N); 
Yl=interp2(posY,posX,Zl,posY0,posX0,'bicubic');

% Extracting features
Eh=Yh-Yl; % pre-processing of the high-resolution image
Yl1=conv2(Yl,[1,0,-1],'same'); % the filter is centered and scaled well for s=3
Yl2=conv2(Yl,[1,0,-1]','same');
Yl3=conv2(Yl,[1,0,-2,0,1]/2,'same');
Yl4=conv2(Yl,[1,0,-2,0,1]'/2,'same');

% Gathering the patches
Ph=zeros(n^2,(N/s-2*dd)^2); 
Ptilde_l=zeros(4*n^2,(N/s-2*dd)^2); 
n2=(n-1)/2; 
counter=1;
for k1=s*dd+1:s:N-s*dd
    for k2=s*dd+1:s:N-s*dd
        Ph(:,counter)=reshape(Eh(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);
        Ptilde_l(:,counter)=[reshape(Yl1(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]); ...
                                    reshape(Yl2(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl3(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl4(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1])];
        counter=counter+1;
    end;
end;

% Dimentionalily reduction
R=Ptilde_l*Ptilde_l'; 
[B,SS]=eig(R); 
Permute=fliplr(eye(size(R,1))); 
SS=Permute*SS*Permute; % so that the eigenvalues are sorted descending
B=B*Permute; 
energy=cumsum(diag(SS))/sum(diag(SS)); 
% figure(1); clf; plot(energy)
pos=find(energy>0.999,1);
B=B(:,1:pos);
disp(['Effective dimension: ',num2str(pos)]); 
disp(['The relative error is: ',...
        num2str(mean(sum((B*B'*Ptilde_l-Ptilde_l).^2))/...
        mean(sum((Ptilde_l).^2)))]); % showing the relative error
Pl=B'*Ptilde_l; 

% Low-Res. Dictionary Learning 
param.errorFlag=0;
param.K=m; 
param.numIteration=20; 
param.InitializationMethod='DataElements'; 
param.TrueDictionary=0;
param.Method='KSVD';
param.L=L;   
[Al,output]=Chapter_12_TrainDic_Fast(Pl,param);
Q=output.CoefMatrix; 

% High-Resolution Dictionary Leanring
Ah=Ph*Q'*inv(Q*Q');

%============================================
%                  Sanity Check - Interpolating the training image
%============================================

Ph_hat=Ah*Q; 
Yout=Yl*0; 
Weight=Yl*0;
counter=1;
for k1=s*dd+1:s:N-s*dd
    for k2=s*dd+1:s:N-s*dd
        patch=reshape(Ph_hat(:,counter),[n,n]);
        Yout(k1-n2:k1+n2,k2-n2:k2+n2)=...
            Yout(k1-n2:k1+n2,k2-n2:k2+n2)+patch; 
        Weight(k1-n2:k1+n2,k2-n2:k2+n2)=...
            Weight(k1-n2:k1+n2,k2-n2:k2+n2)+1; 
        counter=counter+1;
    end;
end;
Yout=Yout./(Weight+1e-5)+Yl; 
figure(2); clf; imagesc([Yh,Yl,Yout]);
colormap(gray(256)); axis image; axis off; 

ErrOut=mean(mean((Yout(dd*s+1:end-dd*s,dd*s+1:end-dd*s)...
                           -Yh(dd*s+1:end-dd*s,dd*s+1:end-dd*s)).^2));
ErrIn=mean(mean((Yl(dd*s+1:end-dd*s,dd*s+1:end-dd*s)...
                           -Yh(dd*s+1:end-dd*s,dd*s+1:end-dd*s)).^2)); 
disp('Result on the train image: ');
disp([sqrt(ErrIn),sqrt(ErrOut),10*log10(ErrIn/ErrOut)]);

% The only variable to keep:Al, Ah, B, and (dd,L,n,n2,m,s)
% clear  Eh Permute Ph Ph_hat Ptilde_l Pl Q R SS Yout param
% clear output counter patch pos k1 posY0 energy N k2
% clear Weight Yh Yl Yl1 Yl2 Yl3 Yl4 Zl posX posY posX0 

%============================================
%                                Interpolating a test image
%============================================
         
Yh=imread(ImageIN);
N=size(Yh,1);
N=floor(N/s^2)*s^2; 
Yh=im2double(Yh(1:N,1:N)); % so that it dcales down to an integer size
Yh=Yh*255; 

% Creating the low-resolution image
Zl=conv2(Yh,[1 2 1]/4,'same');
Zl=conv2(Zl,[1 2 1]'/4,'same');
Zl=Zl(1:s:end,1:s:end); 

% Upscaling Zl to the original resolution
[posY,posX]=meshgrid(1:s:N,1:s:N); 
[posY0,posX0]=meshgrid(1:1:N,1:1:N); 
Yl=interp2(posY,posX,Zl,posY0,posX0,'bicubic');

% Extracting features
Yl1=conv2(Yl,[1,0,-1],'same'); 
Yl2=conv2(Yl,[1,0,-1]','same');
Yl3=conv2(Yl,[1,0,-2,0,1]/2,'same');
Yl4=conv2(Yl,[1,0,-2,0,1]'/2,'same');

% Gathering the patches
Ptilde_l=zeros(4*n^2,(N/s-2*dd)^2); 
counter=1;
for k1=s*dd+1:s:N-s*dd
    for k2=s*dd+1:s:N-s*dd
        Ptilde_l(:,counter)=[reshape(Yl1(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]); ...
                                    reshape(Yl2(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl3(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1]);...
                                    reshape(Yl4(k1-n2:k1+n2,k2-n2:k2+n2),[n^2,1])];
        counter=counter+1;
    end;
end;

% Dimentionalily reduction
Pl=B'*Ptilde_l; 

% Cleaning up some space in memory
clear posX posY posX0 posY0 Zl Yi4 Yl3 Yl2 Yl1 Ptilde_l

% Sparse coding of the low-res patches
% Q=omp(Al'*Pl,Al'*Al,L); does not run because of memory problems
Q1=omp(Al'*Pl(:,1:60000),Al'*Al,L);
Q2=omp(Al'*Pl(:,60001:120000),Al'*Al,L);
Q3=omp(Al'*Pl(:,120001:end),Al'*Al,L);
Q=[Q1,Q2,Q3];

% Recover the image
Ph_hat=Ah*Q; 
Yout=Yl*0; 
Weight=Yl*0;
counter=1;
for k1=s*dd+1:s:N-s*dd
    for k2=s*dd+1:s:N-s*dd
        patch=reshape(Ph_hat(:,counter),[n,n]);
        Yout(k1-n2:k1+n2,k2-n2:k2+n2)=...
            Yout(k1-n2:k1+n2,k2-n2:k2+n2)+patch; 
        Weight(k1-n2:k1+n2,k2-n2:k2+n2)=...
            Weight(k1-n2:k1+n2,k2-n2:k2+n2)+1; 
        counter=counter+1;
    end;
end;
Yout=Yout./(Weight+1e-5)+Yl; 
Yout=min(max(Yout,0),255);
Yout(:,end)=255; % annoying boundary problem (does not influence the MSE)
figure(1); clf; 
imagesc([Yh,ones(N,5)*255, 2*abs(Yl-Yout); ...
              ones(5,2*N+5)*255; 
              Yl,ones(N,5)*255,Yout])
colormap(gray(256)); axis image; axis off; 

figure(2); clf; 
image([Yh(601:700,451:550),ones(100,5)*255,...
           Yl(601:700,451:550),ones(100,5)*255,...
           Yout(601:700,451:550)]); 
colormap(gray(256)); axis image; axis off;
figure(3); clf; 
image([Yh(301:400,201:300),ones(100,5)*255,...
           Yl(301:400,201:300),ones(100,5)*255,...
           Yout(301:400,201:300)]); 
colormap(gray(256)); axis image; axis off;

ErrOut=mean(mean((Yout(dd*s+1:end-dd*s,dd*s+1:end-dd*s)...
                           -Yh(dd*s+1:end-dd*s,dd*s+1:end-dd*s)).^2));
ErrIn=mean(mean((Yl(dd*s+1:end-dd*s,dd*s+1:end-dd*s)...
                           -Yh(dd*s+1:end-dd*s,dd*s+1:end-dd*s)).^2)); 
disp('Result on the test image: ');
disp([sqrt(ErrIn),sqrt(ErrOut),10*log10(ErrIn/ErrOut)]);

return; 











