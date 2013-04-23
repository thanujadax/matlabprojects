function offEdgeIDs = getUnOrientedEdgeIDs(edges2pixels,orientedScoreSpace3D,...
                lenThresh,orientationScoresMax)
% returns the IDs of edges that have sharp changes of orientation (~180)
% along its pixels

% Inputs:
%   

% get the list of edges with the corresponding set of pixels for each edge
% for each edge of each pixel, get the maxRespOrientation
% get the diff of the maxRespOrientation values
% if it has a diff around 180 degrees, discard that edge

[numEdges,maxPixelsPerEdge] = size(edges2pixels);
edgeLengths = zeros(numEdges,1);
for i = 1:numEdges
    edgeLengths(i) = sum((edges2pixels(i,:))>0) - 1;
end

edgeListIndToExamine = find(edgeLengths<lenThresh);
for i = 1:numel(edgeListIndToExamine)
    % check if the edge contains a unacceptable set of orientations
    % if yes, add it to the list of offEdgeIDs
    
end
