function z = getZ(Dictionary,data,verticalAMat,horizontalAmat,N)
% calculate the normalization factor z for a given 4-neighborhood N(i) and
% for a given input patch

% Inputs:
% Dictionary
% data - the input patch for the current position
% verticalAMat - vertical associations
% horizontalAMat - horizontal associations
% N - vector containing the 4 neighborhood elements

% Output:
% z - normalization factor

f1 = getUnaryPotentialVec(data,Dictionary);    % column vector
f2 = getPairwisePotentialVec(Dictionary,N,verticalAMat,horizontalAMat); % column vector

expPotentials = exp(-f1-f2);

z = sum(expPotentials);

