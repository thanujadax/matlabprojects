function plotAllInterpolationCurves(outputSavePath)

% script to plot different calibration curves for the same volume
% outputSavePath contains all the .mat files that encodes the estimation
% curve. These files are created by doThicknessEstimation().
%% Parameters

% calibrationMethod
% 1 - correlation coefficient across ZY sections, along x axis
% 2 - correlation coefficient across XY sections, along x axis
% 3 - SD of XY per pixel intensity difference
% 4 - c.o.c across ZY sections along Y
% 5 - c.o.c across XY sections, along X
% 6 - c.o.c acroxx XZ sections, along X
% 7 - 
% TODO: methods robust against registration problems

params.xyResolution = 5; % nm
params.maxShift = 15;
params.maxNumImages = 600; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
params.numPairs = 2; % number of section pairs to be used to estimate the thickness of onesection
params.plotOutput = 1;
params.usePrecomputedCurve = 0;
params.pathToPrecomputedCurve = '';

inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes/s108/s108.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s502';

fileStr = 'xcorrMat'; % general string that defines the .mat file

%% Plot all interpolation curves with shaded error bars with different colors
% x
% y - C x N matrix. C is the number of curves. N is the number of samples
% errBar - vector with SD.

% figure; title('Thickness interpolation curves');
% mseb(x,y_mean,y_std);ylim([-50 150])
% legend('Line 1','Line 2','Line 3')
% figure; title('openGL');
% mseb(x,y_mean,y_std,[],1);ylim([-50 150]) 
% legend('Line 1','Line 2','Line 3')

x = 1:params.maxShift;
[y,errBar] = getMeanInterpolationCurves(outputSavePath,fileStr);
lineProps = [];
transparent = 1;

figure; 
H = mseb(x,y,errBar,lineProps,transparent);
legend('1.ZY_x','2.XY_y','4.ZY_y','5.XY_x','6.XZ_x','7.XZ_y');
xlabel('Distance (num pixels)')
ylabel('Coefficient of correlation')
title('Thickness interpolation curves: s502');
%set(gca,'position',[0 0 1 1],'units','normalized')

% same plot without error bars
figure; 
% errBarZero = zeros(size(errBar));
% H2 = mseb(x,y,errBarZero,lineProps,transparent);
plot(y');
legend('1.ZY_x','2.XY_y','4.ZY_y','5.XY_x','6.XZ_x','7.XZ_y');
xlabel('Distance (num pixels)')
ylabel('Coefficient of correlation')
title('Thickness interpolation curves: s502');
%set(gca,'position',[0 0 1 1],'units','normalized')
%% Plot all predicted thicknesses
predictedThickness = getPredictedThicknessesFromTxtFile(outputSavePath);
errBarZ = zeros(size(predictedThickness));
figure;
% H3 = mseb([],predictedThickness,errBarZ,lineProps,transparent);
plot(predictedThickness')
legend('1.ZY_x','2.XY_y','4.ZY_y','5.XY_x','6.XZ_x','7.XZ_y');
xlabel('Section interval index')
ylabel('Estimated thickness (nm)')
title('Thickness estimates: s502')
%set(gca,'position',[0 0 1 1],'units','normalized')
%% Calculate variance of prediction across different methods, from the same section
