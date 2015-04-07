function relZresolution = predictThicknessFromCurve(...
        imageStackFileName,xcorrMat,maxShift)
% Returns the section thickness relative to the xy resolution. Multiply by
% xyResolution to get the actual thickness.
    
inputImageStack = readTiffStackToArray(imageStackFileName);
numImg = size(inputImageStack,3);

relZresolution = zeros(1,numImg-1); % relative to xy pix resolution

for i = 1:numImg-1
   
   image1 = inputImageStack(:,:,i);
   image2 = inputImageStack(:,:,(i+1));
   relZresolution(i) = getRelativeDistance_cc2(image1,image2,mean(xcorrMat,1),maxShift);
   
end