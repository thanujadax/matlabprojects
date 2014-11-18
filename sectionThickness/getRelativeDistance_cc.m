function relativeDistancePix = getRelativeDistance_cc(imageName1,imageName2,curve,maxVal)
% predicts relative section interval in terms of the pixel length in xy
% plane, using correlation coefficient

I1 = double(imread(imageName1));
I2 = double(imread(imageName2));

cc = corr2(I1,I2);

relativeDistancePix = interp1(curve,1:maxVal,cc);
