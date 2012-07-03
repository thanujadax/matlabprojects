%inputfile = '/home/thanuja/matlabprojects/data/mitoData/stem1.tiff';
inputfile = '/home/thanuja/matlabprojects/data/mitoData/schmidhuber1.tif';

I=imread(inputfile,'tiff');


%mask = [-1 -1 -1;-1 8 -1;-1 -1 -1;];

mask = [1 0 0 0 0 0 0 0 0 0 0 0 0 0
     1 1 0 0 0 0 0 0 0 0 0 0 0 0
     1 1 1 0 0 0 0 0 0 0 0 0 0 0
     1 1 1 1 0 0 0 0 0 0 0 0 0 0 
     1 1 1 1 1 0 0 0 0 0 0 0 0 0
     1 1 1 1 1 1 1 0 0 0 0 0 0 0 
     1 1 1 1 1 1 1 1 1 1 0 0 0 0 
     1 1 1 1 1 1 1 1 1 1 1 1 1 1];

figure();
imshow(I);

%J = conv2(I, mask, 'valid');
J = xcorr2(I, mask);



figure();
imshow(J); 