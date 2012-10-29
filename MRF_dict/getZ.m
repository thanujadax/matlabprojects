function z = getZ(Dictionary,data,verticalAMat,horizontalAmat,currentLabels,rowSize,colSize)
% calculate the normalization factor z for a given 4-neighborhood N(i) and
% for a given input patch

% Inputs:
% Dictionary
% data - the input patch for the current position
% verticalAMat - vertical associations
% horizontalAMat - horizontal associations
% rowSize - number of patches per row in the image
% colSize - number of patches per column in the image

% Output:
% z - normalization factor

f1 = getUnaryPotentialVec(data,Dictionary);    % column vector
f2 = getPairwisePotentialVec(Dictionary,currentLabels,verticalAMat,horizontalAMat,rowSize,colSize); % column vector

expPotentials = exp(-f1-f2);

z = sum(expPotentials);

