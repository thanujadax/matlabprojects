function peaks3D = houghBarPeaks(houghSpace3D,orientations,thresholdFraction...
                    ,slidingDist,barLength,barWidth)

% TODO: parallelize                
                
% output:
%   peaks3D - [row col orientation] 3D array containing the votes for the
%   detected peaks or zero otherwise.

% inputs:
%   houghSpace3D - contains the vote for each pixel at each orientation as
%   a 3D array [row col orientation]
%   orientations - vector of the orientations e.g. [0 45 90 135]
%   thresholdFraction - what fraction of max(HoughVote) should be used for
%   thresholding the peaks
%   slidingDist - spacing between pixels for voting 
%   barLength -
%   barWidth - 

% 1. set the threshold for peak detection separately per each orientation's max vote
% 2. extract all possible peaks for each orientation (thresholding)
% 3. non-max suppression in the defined neghborhoods
% 4. 

peaks3D = zeros(size(houghSpace3D));
[numRows numCols numOrientations] = size(houghSpace3D);

for i = 1:numOrientations
    % get max vote for this orientation
    maxVote = max(max(houghSpace3D(:,:,i)));
    thresh = maxVote*thresholdFraction;
    [r,c,vote] = find(houghSpace3D(:,:,i)>thresh);
    voteInd = sub2ind([numRows numCols],r,c);
    voteMat = zeros(numRows,numCols);
    voteMat(voteInd) = vote;
    
    peaks3D(:,:,i) = NMS_bars(voteMat,orientations(i),barLength,barWidth);
    
end