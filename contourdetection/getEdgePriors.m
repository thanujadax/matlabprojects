function [edgePriors,edgeOrientationIndex] = getEdgePriors(orientedScoreSpace3D,edges2pixels)
% computes N-by-nOrientation matrix for each edge given in edges2pixels
% returns the max response for each edge
% for each edge, get the edge pixels. then get the orientation response to
% all the pixels. For each orientation, calculate the mean response. choose
% the orientation with the max avg response and the corresponding response
% as the edgePrior.
% Also return the edgeOrientationIndex as the 2nd output. To get the actual
% edgeOrientation, multiply the orientationInd by the orientationStepSize.

% Inputs:
%   orientedScoreSpace3D - m-by-n-by-nOrientation matrix for the
%   orientation response for each pixel in the image of size m-by-n
%   edges2pixels - contains the pixel inds for each edge. first col:
%   edgeIDs

% Output:
%   edgePriors - N-by-1 array of max responses
[~,nce] = size(edges2pixels);  % first column is the edgeID
edgepixels = edges2pixels(:,2:nce);
[~,~,nOrientations] = size(orientedScoreSpace3D); 
numEdges = size(edgepixels,1);
edgePriors_all = zeros(numEdges,nOrientations);
edgePriors = zeros(numEdges,1);
edgeOrientationIndex = zeros(numEdges,1);
for i=1:numEdges
    % for each edge, get the pixel indices
    edgePixelInds = edgepixels(i,:);              % list indices
    edgePixelInds = edgePixelInds(edgePixelInds>0); % indices of edge pixels wrt image
    if(~isempty(edgePixelInds))
        % for each orientation take the average over all pixels
        for j=1:nOrientations
            % get the total response for all edge pixels for this dimension
            orientationResp_j = orientedScoreSpace3D(:,:,j);
            meanEdgePixResp_j = mean(orientationResp_j(edgePixelInds));
            edgePriors_all(i,j) = meanEdgePixResp_j;
        end
        [edgePriors(i),edgeOrientationIndex(i)] = max(edgePriors_all(i,:));  
    else
        edgePriors(i) = 0;  % self edge
    end
end