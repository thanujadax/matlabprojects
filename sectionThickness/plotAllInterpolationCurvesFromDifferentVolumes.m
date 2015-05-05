function plotAllInterpolationCurvesFromDifferentVolumes(outputSavePath)

% script to plot different calibration curves for different volumes
% outputSavePath contains all the .mat files that encodes the estimation
% curve. These files are created by doThicknessEstimation().
% usually, the outputSavePath contains the calibration curves for different
% volumes using the same calibration method.
%% Parameters

% calibrationMethod
% % 1 - c.o.c across ZY sections, along x axis
% % 2 - c.o.c across XY sections, along Y axis
% % 3 - SD of XY per pixel intensity difference
% % 4 - c.o.c across ZY sections along Y
% % 5 - c.o.c across XY sections, along X
% % 6 - c.o.c acroxx XZ sections, along X
% % 7 - c.o.c acroxx XZ sections, along Y 
% TODO: methods robust against registration problems

calibrationMethod = 1;

params.xyResolution = 5; % nm
params.maxShift = 15;
params.maxNumImages = 100; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
params.numPairs = 1; % number of section pairs to be used to estimate the thickness of onesection
params.plotOutput = 1;
params.usePrecomputedCurve = 0;
params.pathToPrecomputedCurve = '';

outputSavePath = '/home/thanuja/projects/tests/thickness/batchEstimation/20150504/calibration01';

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

% titleStr = 'Thickness interpolation curves';
titleStr = sprintf('Thickness interpolation curves: calibration %d',calibrationMethod);
figure; title(titleStr);
H = mseb(x,y,errBar,lineProps,transparent);
legend('s108','s202','s410','s502','s603','s704','s801','s901','s909');
xlabel('Distance (num pixels)')
ylabel('Coefficient of correlation')

% same plot without error bars
figure; 
% errBarZero = zeros(size(errBar));
% H2 = mseb(x,y,errBarZero,lineProps,transparent);
plot(y');
title(titleStr);
legend('s108','s202','s410','s502','s603','s704','s801','s901','s909');
xlabel('Distance (num pixels)')
ylabel('Coefficient of correlation')

