function z = getZ(Dictionary,inputData,verticalAMat,horizontalAMat,currentLabels,rowSize,colSize)
% calculate the normalization factor z for a given 4-neighborhood N(i) and
% for a given input patch

% Inputs:
% Dictionary
% data - the input patch for the current position. if data = 0, no input is
% assumed. i.e there will not be any unary potential calculated.
% verticalAMat - vertical associations
% horizontalAMat - horizontal associations
% rowSize - number of patches per row in the image
% colSize - number of patches per column in the image

% Output:
% z - normalization factor

if(inputData~=0)
    f1 = zeros(1,size(Dictionary,1));
else
    f1 = getUnaryPotentialVec(inputData,Dictionary);    % column vector
end

f2 = getPairwisePotentialVec(currentLabels,verticalAMat,horizontalAMat,rowSize,colSize); % column vector

expPotentials = exp(-f1-f2);

z = sum(expPotentials);

