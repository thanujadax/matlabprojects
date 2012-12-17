function patch = drawLineInMat(patch,x,y)
% Inputs:
%   patch - input matrix (which will be edited)
%   x - vector of x coordinates of the points
%   y - vector of y coordinates of the points

matSize = size(patch); 

% convert the vector of subscripts into a vector of indices
x = int16(x);
y = int16(y);
indices = sub2ind(matSize,y,x);

% set the corresponding elements to 1
patch(indices) = 1;
