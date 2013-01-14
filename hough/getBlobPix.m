function blobPix = getBlobPix(peaks3D,orientations,blobThresh)

% Inputs:
%   peaks3D: peaks of the 3D hough space [r c orientation]
%   orientations: vector of the orientations used in degrees
%   blobThresh: threshold fraction of max(peak3D) for blob detection

[numRows numCols numOrientations] = size(peaks3D);

% we have the peaks for each orientation

% get the internal blob points i.e. the points having a high vote in all
% orientations
blobPix = zeros(numRows,numCols,numOrientations);

parfor i = 1:numOrientations
    blobMat = zeros(numRows,numCols);
    voteMat = peaks3D(:,:,i);
    maxVote = max(max(voteMat));
    thresh = maxVote * blobThresh;
    internalBlobInd = find(voteMat>thresh);
    blobMat(internalBlobInd) = 1;
    blobPix(:,:,i) = blobMat;    
end

blobVotes = sum(blobPix,3);
% starting from the internal blob points, go towards the edge to get the
% boundary points of the blob