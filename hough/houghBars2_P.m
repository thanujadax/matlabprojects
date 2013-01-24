function houghSpace3D = houghBars2_P(img,barLength,barWidth,orientations,slidingDist)

% inputs: 
%   img - contains the pixels that should be used for calculating the hough
%   space

% output:
% houghSpace3D

% each point in concern, votes for all the bars that it is contributing, at once
[numRows numCols] = size(img);
numOrientations = size(orientations,2);

% init
houghSpace3D = zeros(numRows,numCols,numOrientations);

margin = ceil(max(barLength,barWidth)+1/2);   % allowing a margin in the image

% set the pixels outside this margin to zero
% img(1:margin,:) = 0;
% img((numRows-margin):numRows,:) = 0;
% img(:,1:margin) = 0;
% img(:,(numCols-margin):numCols) = 0;

% progressbar('Calculating 3D Hough space'); % Create figure and set starting time
% get the indices of the pixels that should be voted with
inds = find(img>0.4);
totPoints = numel(inds);

parfor i=1:numOrientations
    inds_i = inds;
    img_i = img;
    for j=1:totPoints
        ind = inds_i(j);
        pixVal = img_i(ind);
        % get the indices of the bar centers that can include this pixel
        [r c] = ind2sub([numRows numCols],ind);
        barCentInd = getBarPixInd(r,c,orientations(i),barLength,barWidth,numRows,numCols);
        houghVote2D_i = zeros(numRows,numCols);
        houghVote2D_i(barCentInd) = pixVal;
        houghSpace3D(:,:,i) = houghSpace3D(:,:,i) + houghVote2D_i;
    end
  % progressbar(j/totPoints); % update progress bar
end