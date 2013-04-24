function onEdgeListIDs = getBackboneEdgeIDs(edgepixels,edgePriors,...
                lenThreshBB,priorThreshFracBB)
% returns the IDs of edges that have a prior response higher than the
% threshold, so that these edges can be considered as the backbone of the
% resulting segmentation.

% Inputs:
%   

% get the list of edges higher than the min length (lenThreshBB)
% out of those edges, pick the edges that have an edgePrior above threshold
% add them to the onEdgeListIDs

[numEdges,maxPixelsPerEdge] = size(edgepixels);
edgeLengths = zeros(numEdges,1);
for i = 1:numEdges
    edgeLengths(i) = sum((edgepixels(i,:))>0);
end

edgeListIndToExamine = find(edgeLengths>lenThreshBB);
edgePriors_examineList = edgePriors(edgeListIndToExamine);

onEdgeListIDs = edgeListIndToExamine(edgePriors_examineList>priorThreshFracBB);