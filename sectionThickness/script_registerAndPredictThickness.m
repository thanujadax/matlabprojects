% script register and predict thickness

inputImagePath = '/home/thanuja/projects/inputData/trainingHalf/raw';
imgFileType = 'png';
patchSizeX = 128;
patchSizeY = 128;
maxNumPatches = 20;
overlap = 20;
maxThicknessPix = 25;

[t_median,t_mean,t_var, thicknessCurve,tCurve_std] ...
    = registerAndPredictThickness(...
    inputImagePath,imgFileType,patchSizeX,patchSizeY,maxNumPatches,overlap,...
    maxThicknessPix);