function [internalx,internaly] = getInternelPixelsFromBoundary(boundaryPixels,sizeR,sizeC)

% Inputs:
%   boundaryPixels - vector of the pixelinds of the boundary pixels
%   sizeR - num rows of original image
%   sizeC - num cols of original image

[y,x] = ind2sub([sizeR sizeC],boundaryPixels);

multipleHighest = 0;        % set flag if multiple points are found at the top
multipleLowest = 0;         % set flag if multiple points are found at the bottom

[~,highestpt_listInd] = max(y);
[~,lowestpt_listInd] = min(y);

% if multiple maxima/minima found, pick one of them
if(numel(highestpt_listInd)>1)
    highestpt_listInd = highestpt_listInd(1);
    multipleHighest =1;
end
if(numel(lowestpt_listInd)>1)
    lowestpt_listInd = lowestpt_listInd(1);
    multipleLowest = 1;
end

highestpt.y = y(highestpt_listInd);
highestpt.x = x(highestpt_listInd);

lowestpt.y = y(lowestpt_listInd);
lowestpt.x = x(lowestpt_listInd);

numHorizontalLines = highestpt.y - lowestpt.y - 1; % excluding the boundaries

internalx = [];
internaly = [];

for i=(lowestpt.y+1):(highestpt.y-1)
    % i is the horizontal line number i.e. the y coordinate
    % get the left most and right most points for this horizontal line
    pixListInd_i = find(y==i);
    % get all x coordinates with this y values
    x_for_y = x(pixListInd_i);
    x_left = min(x_for_y);          % left most x
    x_right = max(x_for_y);         % right most x
    numRowPixels = x_right - x_left - 1; % excluding the boundary pixels
    j=(x_left+1):(x_right-1);        % all internal x coordinates for y=i
    k = ones(numRowPixels,1) .* i;  % the y coordinate vector
    
    % append to the output column vectors
    internalx = [internalx; j'];
    internaly = [internaly; k];
end