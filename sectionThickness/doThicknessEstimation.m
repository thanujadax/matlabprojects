function thicknessEstimates = doThicknessEstimation(...
    calibrationMethod,inputImageStackFileName,outputSavePath,params)

% Performs section thickness estimation using representative
% curves to determine the distance between two (adjacent) sections
% do thickness estimation based on one of the following methods to
% calibrate the similarity curve
%
% output:
%   thicknessEstimates - column vector of thickness estimates = distance
%   between 2 adjacent sections. It will have multiple columns if
%   param.numPairs > 1. 2nd column will be the estimate using the i and
%   i+2nd image. etc.

%% Parameters

% calibrationMethod = 5;
% 
% % 1 - c.o.c across ZY sections, along x axis
% % 2 - c.o.c across XY sections, along Y axis
% % 3 - SD of XY per pixel intensity difference
% % 4 - c.o.c across ZY sections along Y
% % 5 - c.o.c across XY sections, along X
% % 6 - c.o.c acroxx XZ sections, along X
% % 7 - c.o.c acroxx XZ sections, along Y
% % TODO: methods robust against registration problems
% 
% params.predict = 1; % set to 0 if only the interpolation curve is required.
% params.xyResolution = 5; % nm
% params.maxShift = 15;
% params.maxNumImages = 60; % number of sections to initiate calibration.
%                 % the calibration curve is the mean value obtained by all
%                 % these initiations
% params.numPairs = 1; % number of section pairs to be used to estimate the thickness of onesection
% params.plotOutput = 1;
% params.usePrecomputedCurve = 0;
% params.pathToPrecomputedCurve = '';
% 
% inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes/s108/s108.tif';
% outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s108';

%% 

if(calibrationMethod == 1)
    disp('Calculating c.o.c decay curve using ZY stack, along X ...')
    calibrationString = 'c.o.c ZY along X';
    calibrationFigureFileString = '01_cocZY_X';
    xcorrMat = getXcorrZYstack(inputImageStackFileName,params.maxShift,params.maxNumImages);
    disp('done!')
    % each column of xcorrMat corresponds to a sequence of shifted frames of
    % the zy plane along the x axis.
    % each row corresponds to one starting image (zy section) of the stack

elseif(calibrationMethod == 2)
    disp('Calculating c.o.c decay curve using XY images stack, along Y ...')
    calibrationString = 'c.o.c XY along Y';
    calibrationFigureFileString = '02_cocXY_Y';
    xcorrMat = getXcorrXYstack(inputImageStackFileName,params.maxShift,params.maxNumImages);
    disp('done!')    
    
elseif(calibrationMethod == 3)
    calibrationString = 'SD of pixel intensity XY along X';
    disp('Calculating SD of intensity deviation curve using shifted XY sections stack, along X ...')
    xcorrMat = getIntensityDeviationXYstack(inputImageStackFileName,params.maxShift,params.maxNumImages);
    disp('done!')
    calibrationFigureFileString = '03_sdpiXY_X';
    
elseif(calibrationMethod == 4)
    calibrationString = 'c.o.c ZY along Y';
    disp('Calculating c.o.c decay curve using ZY stack, along Y ...')
    xcorrMat = getXcorrZYstackY(inputImageStackFileName,params.maxShift,params.maxNumImages);
    disp('curve estimation done')
    calibrationFigureFileString = '04_cocZY_Y';
    
elseif(calibrationMethod == 5)
    calibrationString = 'c.o.c XY along X';
    disp('Calculating c.o.c decay curve using XY images stack, along X ...')
    xcorrMat = getXcorrXYstackX(inputImageStackFileName,params.maxShift,params.maxNumImages);
    disp('curve estimation done!')    
    calibrationFigureFileString = '05_cocXY_X';
    
elseif(calibrationMethod == 6)
    calibrationString = 'c.o.c XZ along X';
    disp('Calculating c.o.c decay curve using XZ images stack, along X ...')
    xcorrMat = getXcorrXZstackX(inputImageStackFileName,params.maxShift,params.maxNumImages);
    disp('curve estimation done!')    
    calibrationFigureFileString = '06_cocXZ_X';
    
elseif(calibrationMethod == 7)
    calibrationString = 'c.o.c XZ along Y';
    disp('Calculating c.o.c decay curve using XZ images stack, along Y ...')
    xcorrMat = getXcorrXZstackY(inputImageStackFileName,params.maxShift,params.maxNumImages);
    disp('curve estimation done!')
    calibrationFigureFileString = '07_cocXZ_Y';
    
else
    error('Unrecognized calibration method specified. Check calibrationMethod')
end


%% plot calibration curve

tokenizedFName = strsplit(inputImageStackFileName,filesep);
nameOfStack = strtok(tokenizedFName(end),'.');
nameOfStack = nameOfStack{1};

if(params.plotOutput)
    titleStr = sprintf('Similarity curve (%d): %s. Vol %s',calibrationMethod,calibrationString,nameOfStack);
    title(titleStr)
    xlabelStr = 'Shifted pixels';
    ylabelStr = 'Coefficient of Correlation';
    transparent = 0;
    shadedErrorBar((1:params.maxShift),mean(xcorrMat,1),std(xcorrMat),'g',transparent,...
        titleStr,xlabelStr,ylabelStr);
    % save calibration curve figure
    % calibrationFigureFileName = strcat(calibrationFigureFileString,'.png')
    calibrationFigureFileName = fullfile(outputSavePath,calibrationFigureFileString);
    print(calibrationFigureFileName,'-dpng');
end

%% save calibration curve .mat
disp('Saving xcorrMat')

matName = sprintf('xcorrMat_%s_%02d.mat',nameOfStack,calibrationMethod);
xcorrFileName = fullfile(outputSavePath,matName);
save(xcorrFileName,'xcorrMat')
disp('done')

%% Predict thickness
if(params.predict)
% predict section thickness for the data set
% relZresolution = predictThicknessFromCurve(...
%         inputImageStackFileName,xcorrMat,params.maxShift,calibrationMethod);
    
relZresolution = predictThicknessFromCurveFromMultiplePairs(...
        inputImageStackFileName,xcorrMat,params.maxShift,calibrationMethod,params.numPairs);
% each row contains one set of estimates for all the sections. First row is
% from the images i and i+1, the second row is from images i and i+2 etc    
thicknessEstimates = relZresolution .* params.xyResolution;
%% save predicted thickness to text file
thicknessFileName = sprintf('_cm%02d_thickness.txt',calibrationMethod);

thicknessFileName = strcat(nameOfStack,thicknessFileName);
thicknessFileName = fullfile(outputSavePath,thicknessFileName);
thicknessEstimates = thicknessEstimates';
save(thicknessFileName,'thicknessEstimates','-ASCII');
%% plot output if required
if(params.plotOutput)
    if(size(thicknessEstimates,2)==1)
        figure;plot(thicknessEstimates);
    elseif(size(thicknessEstimates,2)==2)
        figure;plot(thicknessEstimates(:,1),'r');
        hold on
        plot((thicknessEstimates(:,2).* 0.5),'g');
        hold off
    end
    titleStr = sprintf('Section thickness estimates. Vol: %s. Cal: %d',nameOfStack,calibrationMethod);
    title(titleStr);
    xlabel('Section index')
    ylabel('Estimated thickness (nm)')
end
end
