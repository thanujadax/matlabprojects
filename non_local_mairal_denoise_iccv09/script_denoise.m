images={'house.png','peppers256.png','Cameraman256.png','lena.png','barbara.png','boat.png','hill.png','couple.png','man.png','fingerprint.png','bridge.png','flintstones.png'}; 
tabsigma=[5 10 15 20 25 50 100];
format compact;

if (exist('res.mat'))
   load('res.mat');
else
   res=zeros(length(tabsigma),length(images));
end

for ii=1:length(images);
   I=double(imread(sprintf('data/%s',images{ii})))/255;
   if (size(I,3) > 1)
      I=I(:,:,1);
   end
   Io=I;
   for jj=1:length(tabsigma)
      if (res(jj,ii) == 0)
         randn('seed',0);
         sigma=tabsigma(jj)/255;
         I=Io+(sigma)*randn(size(I));
         Iout=denoise(I,sigma,Io);
         res(jj,ii)=psnr(Io,Iout);
         res
         save('res.mat','res')
      end
   end
end
