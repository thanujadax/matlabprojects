function plotMembraneProfiles()

inputImageFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/smallSet1/1.tif';

inputImageMat = double(imread(inputImageFileName));
inputImageMat = inputImageMat./255;

% rescale from -1 to +1
inputImageMat_1_1 = rescaleImg1_1(inputImageMat);

%% membrane 1
rstart =381;
rstop = 398;
col = 252;
vProf = inputImageMat_1_1(rstart:rstop,col);
x = 1:length(vProf);
figure;plot(x,vProf);title('membrane intensity profile');xlabel('pixel id');ylabel('scaled intensity');
%% membrane 2
rstart = 133;
rstop = 146;
col = 309;
vProf = inputImageMat_1_1(rstart:rstop,col);
x = 1:length(vProf);
figure;plot(x,vProf);title('membrane intensity profile');xlabel('pixel id');ylabel('scaled intensity');
