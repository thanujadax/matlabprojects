function [edgeOrientations, edgeResponses] = getEdgeOrientations_grid...
            (edges2nodes,sizeR,sizeC,eps_orientation,OFR,...
            edges2pixels,gridResolution,barHalfWidth,ofr_stepSize)

% Outputs:
%   edgeOrientations (theta): col1 gives default orientation in the order of nodes
%   given in edges2nodes. col2 gives the complementary orientation
%   edgeResponses: edge response in rectangular region in the proximity
%   of the mid point of the edge for the given edge orientation +- epsilon,
%   for the default edgeOrientation (col1 of edgeOrientations). Sign of the
%   edge response suggest whether to use col1(+) or col2(-) in the case of
%   undirected graph version of ILP.

% Inputs:
%   eps_orientation: tolerance on edge orientation for accumulating edge
%   responses
%   OFR: 3D matrix containing edge responses in the order of
%   ofrOrientations
%   OFR matrix:
%   edges2pixels: first col is edgeID. rest of the row gives pixInds
%   gridResolution: dimension (length) of one side of a grid cell
%   barHalfWidth: padding to the edge on one side to calculate edge
%   coverage for edge response calculation
%   ofr_stepSize: difference between adjacent orientation dimensions of ofr in degrees

numEdges = size(edges2nodes,1);

edgeOrientations = zeros(numEdges,2);
edgeResponses = zeros(numEdges,1);

for i=1:numEdges
    % calculation of default edge orientation for each edge, in the
    % direction of node1->node2 as given by edges2nodes
    nodesForEdge = edges2nodes(i,:);
    [theta,theta_comp] = calculateEdgeOrientaion(nodesForEdge,sizeR,sizeC);   
    edgeOrientations(i,1) = theta;
    edgeOrientations(i,2) = theta_comp;
    
    edgesPixels = edges2pixels(i,:);
    edgePixels(1) = []; % first element is the edgeID
    edgePixels = edgePixels(edgePixels>0);
    edgeLength = gridResolution;
    
    edgeOrientation = theta;
    
    pixIndOfr = getEdgeCoverage_rect(edgePixels,edgeLength,...
            sizeR,sizeC,edgeOrientation,barHalfWidth);
    % calculating edge response for default orientation
    edgeResponses(i) = calculateEdgeResp(pixIndOfr,OFR,...
                edgeOrientation,eps_orientation,ofr_stepSize);
    
end



function [theta,theta_comp] = calculateEdgeOrientaion(nodesForEdge,sizeR,sizeC)
% calculation of default edge orientation for each edge, in the
% direction of node1->node2 as given by nodesForEdge = [node1,node2]
[rNodes,cNodes] = ind2sub([sizeR sizeC],nodesForEdge);
y = rNodes(2) - rNodes(1);
x = cNodes(2) - cNodes(1);
theta = atan2d(y,x);        % default orientation
if(theta<0)
    theta = theta + 360;
end

if(theta<180)
    theta_comp = theta + 180;
else
    theta_comp = theta - 180;
end