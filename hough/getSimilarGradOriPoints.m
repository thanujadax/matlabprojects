function gradPointInd = getSimilarGradOriPoints(gradOriMap,orientation)

% parameters
error = 6; % degrees

[numRows numCols] = size(gradOriMap);
orientation_mat = ones(numRows,numCols)*orientation;

diff_mat = sqrt((gradOriMap - orientation_mat)^2);

diff

gradPointInd = find(diff_mat<error);

