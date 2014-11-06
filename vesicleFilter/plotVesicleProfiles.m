% vesicle profiles

inputImageFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/smallSet1/1.tif';

inputImageMat = double(imread(inputImageFileName));
inputImageMat = inputImageMat./255;
% gaussian filtering
sigma = 0.6;
maskSize = 3;
inputImageMat = gaussianFilter(inputImageMat,sigma,maskSize);
inputImageMat_1_1 = rescaleImg1_1(inputImageMat);
figure; imshow(inputImageMat_1_1); title('input image')

%% vesicle 1
rstart =122;
rstop = 147;
col = 255;
vProf = inputImageMat_1_1(rstart:rstop,col);
x = 1:length(vProf);
figure;plot(x,vProf);title('vesicle intensity profile');xlabel('pixel id');ylabel('scaled intensity');

%% vesicle 2
rstart =209;
rstop = 229;
col = 153;
vProf = inputImageMat_1_1(rstart:rstop,col);
x = 1:length(vProf);
figure;plot(x,vProf);title('vesicle intensity profile');xlabel('pixel id');ylabel('scaled intensity');

%% vesicle 3
rstart =194;
rstop = 214;
col = 144;
vProf = inputImageMat_1_1(rstart:rstop,col);
x = 1:length(vProf);
figure;plot(x,vProf);title('vesicle intensity profile');xlabel('pixel id');ylabel('scaled intensity');

%% vesicle 4
rstart =155;
rstop = 175;
col = 304;
vProf = inputImageMat_1_1(rstart:rstop,col);
x = 1:length(vProf);
figure;plot(x,vProf);title('vesicle intensity profile');xlabel('pixel id');ylabel('scaled intensity');