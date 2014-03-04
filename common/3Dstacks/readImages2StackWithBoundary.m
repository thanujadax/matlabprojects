function imageStack3D = readImages2StackWithBoundary...
                (pathForInputImages,fileNameString,gridResY,gridResX)

inputSections = dir(strcat(pathForInputImages,fileNameString));
numInputSections = length(inputSections);
% load images to 3D matrix
% read the first file to read the dimensions
inputImage1_FilePath = fullfile(pathForInputImages,inputSections(1).name);
im1 = uint8(imread(inputImage1_FilePath)); % uint8,single,double are the available options
[numR,numC,sizeZ] = size(im1);
clear im1
% Adjusting image dimensions to include boundary (membrane) cells
numR = numR + gridResY*2;
numC = numC + gridResX*2;
numZ = numInputSections + 2; % Thickness of slices not yet considered
if(sizeZ==1)
    imageStack3D = uint8(ones(numR,numC,numZ));
else
    imageStack3D = uint16(ones(numR,numC,numZ));
end
rStart = gridResY + 1;
rStop = numR - gridResY;
cStart = gridResX + 1;
cStop = numC - gridResX;
for i=1:numInputSections
    % section 1 is already filled with 1s (membrane)
    inputImageFilePath = fullfile(pathForInputImages,inputSections(i).name);
    inputImage = uint8(imread(inputImageFilePath));
    % convert rgb to indexed image if required
    if(sizeZ>1)
        % replace colors with numbers to make RGB_3Dmat->RGB_2Dmat
        [inputImage,numLabels] = getLabelIndexImg_noReorder(inputImage);
    end
    imageStack3D(rStart:rStop,cStart:cStop,(i+1)) = ...
                inputImage; 
    % uint8,single,double are the available options
end