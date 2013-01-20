function gradPointInd = getSimilarGradOriPoints(gradOriMap,orientation)

% inputs
%   margin - the border of the image to leave out

% parameters
error = 6; % degrees

[numRows numCols] = size(gradOriMap);
orientation_mat = ones(numRows,numCols)*orientation;

diff_mat = sqrt((gradOriMap - orientation_mat)^2);

gradPointInd = find(diff_mat<error);

