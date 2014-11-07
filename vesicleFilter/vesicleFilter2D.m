function vesicleFilter2D()

% Vesicle filter 2D

%% Parameters
r1 = 3;
r2 = 6;
r3 = 9;
outputType = 'png';
%% File paths
inputImageFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/smallSet1/1.tif';
outputImageFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/vesicleFilter2d/1.png';
%% Read input
message = sprintf('Input image file name = %s', inputImageFileName);
disp(message)
inputImageMat = double(imread(inputImageFileName));
inputImageMat = inputImageMat./255;
inputImageMat_1_1 = rescaleImg1_1(inputImageMat);
figure; imshow(inputImageMat_1_1); title('input image')


%% Convolution
disp('generate vesicle template')
vesicleTemplate = getVesicleElement2D(r1,r2,r3);
disp('convoluting...')
convolutionResult = convn(inputImageMat_1_1,vesicleTemplate);
figure;imagesc(convolutionResult);title('convolution result')
%% Visualization
message = sprintf('writing output to %s',outputImageFileName);
disp(message);
imwrite(convolutionResult,outputImageFileName,outputType);