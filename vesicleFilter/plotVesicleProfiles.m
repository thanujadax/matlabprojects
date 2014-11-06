% vesicle profiles

inputImageFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/smallSet1/1.tif';

inputImageMat = double(imread(inputImageFileName));
inputImageMat = inputImageMat./255;
inputImageMat_1_1 = rescaleImg1_1(inputImageMat);
figure; imshow(inputImageMat_1_1); title('input image')

% vesicle 1
rstart =122;
rstop = 147;
col = 255;
vProf = inputImageMat_1_1(rstart:rstop,col);
x = 1:len(vProf);
figure;plot(x,vProf)