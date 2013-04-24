function offEdgeListIDs = getUnOrientedEdgeIDs(edges2pixels,orientedScoreSpace3D,...
                lenThresh,orientationScoresMax)
% returns the IDs of edges that have sharp changes of orientation (~180)
% along its pixels

% Inputs:
%   

% get the list of edges with the corresponding set of pixels for each edge
% for each edge of each pixel, get the maxRespOrientation
% get the diff of the maxRespOrientation values
% if it has a diff around 180 degrees, discard that edge

%% parameters
upThresh = 200;
downThresh = 160;
%%
offEdgeListIDs = 0;
[numEdges,maxPixelsPerEdge] = size(edges2pixels);
edgeLengths = zeros(numEdges,1);
for i = 1:numEdges
    edgeLengths(i) = sum((edges2pixels(i,:))>0) - 1;
end

edgeListIndToExamine = find(edgeLengths<lenThresh);
k = 0;
for i = 1:numel(edgeListIndToExamine)
    clear edgePixels_i;
    edgeListID = edgeListIndToExamine(i); % row number for the edge in edges2pixels
    edgePixels_i = edges2pixels(edgeListID,2:maxPixelsPerEdge); % set of pixels for 
                % this edge, padded with zero entries
    edgePixels_i = edgePixels_i(edgePixels_i>0); % remove the zero padding
    pixelOrientations = orientationScoresMax(edgePixels_i).*360; % orientation values for the pixels
    % check if the edge contains a unacceptable set of orientations
    % if yes, add it to the list of offEdgeIDs    
    pxOri_diff = diff(pixelOrientations);
    clear misOrientations;
    misOrientations = intersect(pxOri_diff(pxOri_diff>downThresh),pxOri_diff(pxOri_diff<upThresh));
    % misOrientations should contain an element if this edge is
    % misoriented.
    misOriented = numel(misOrientations);
    if(misOriented>0)
        % edge_i is misoriented
        % add to the off list
        k = k + 1;
        offEdgeListIDs(k) = edgeListID;
    end
end
