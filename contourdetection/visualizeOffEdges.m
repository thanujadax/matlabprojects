function imgOffEdges = visualizeOffEdges(offEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC)
% Inputs
%   offEdgeListIDs - list indices of the edges to be turned off

% Output
%   imgOffEdges - shows which edges are turned off

imgOffEdges = zeros(sizeR,sizeC,3);
nodeMat = zeros(sizeR,sizeC);
edgeMat = zeros(sizeR,sizeC);
offEdgeMat = zeros(sizeR,sizeC);

% select node pixels
nodeMat(nodeInds) = 1;
% get all edge pixels
allEdgePixels = edgepixels(edgepixels>0);

% select all edge pixels
edgeMat(allEdgePixels) = 1;
% select off edge pixels
offEdgePixelMat = edgepixels(offEdgeListIDs,:);
offEdgePixelVec = offEdgePixelMat(offEdgePixelMat>0);
offEdgeMat(offEdgePixelVec) = 1;

% visualization in RGB 
imgOffEdges(:,:,3) = edgeMat;       % blue - all edges
imgOffEdges(:,:,1) = offEdgeMat;    % red - off edges
imgOffEdges(:,:,2) = nodeMat;       % green - nodes
