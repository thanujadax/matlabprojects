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
params.maxNumImages = 60; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
params.numPairs = 2; % number of section pairs to be used to estimate the thickness of onesection
params.plotOutput = 1;
params.usePrecomputedCurve = 0;
params.pathToPrecomputedCurve = '';

inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes/s108.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s108';

fileStr = 'xcorrMat'; % general string that defines the .mat file

%% Plot all curves with shaded error bars with different colors
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

figure; title('Thickness interpolation curves');
H = mseb(x,y,errBar,lineProps,transparent);
legend('1.ZY_x','2.XY_y','4.ZY_y','5.XY_x','6.XZ_x','XZ_y');
