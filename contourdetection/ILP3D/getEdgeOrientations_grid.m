function [edgeOrientations, edgeResponses] = getEdgeOrientations_grid...
            (edges2nodes,sizeR,sizeC,eps_orientation,OFR,ofrOrientations,...
            edges2pixels)

% Outputs:
%   edgeOrientations (theta): col1 gives default orientation in the order of nodes
%   given in edges2nodes. col2 gives the complementary orientation
%   edgeResponses: edge response around a circular region in the proximity
%   of the mid point of the edge for the given edge orientation +- epsilon,
%   for the default edgeOrientation (col1 of edgeOrientations)

% Inputs:
%   eps_orientation: tolerance on edge orientation for accumulating edge
%   responses
%   OFR: 3D matrix containing edge responses in the order of
%   ofrOrientations
%   ofrOrientations: vector containing the angles in the order corresponding to
%   OFR matrix
%   edges2pixels: first col is edgeID. rest of the row gives pixInds

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
    
    % calculating edge response for default orientation
    edgeResp_i = calculateEdgeResp_circular();
    
end

function edgeResp = calculateEdgeResp_circular...
            (edgeID,OFR,ofrOrientations,edges2pixels)
%   calculates edge response around a circular region in the proximity
%   of the mid point of the edge for the given edge orientation +- epsilon,
%   for the default edgeOrientation (col1 of edgeOrientations)

% calculate the relevant pix inds, in the proximity of this edge
edgePixels_i = edges2pixels(edgeID,:);
edgePixels_i(1) = []; % first element is the edgeID
pixIndOfr = getEdgeProximity_circular(edgePixels_i,sizeR,sizeC);
% get the edge response from OFR
edgeResp = calculateWeightedEdgeResp(pixIndOfr,OFR);


function [theta,theta_comp] = calculateEdgeOrientaion(nodesForEdge,sizeR,sizeC)
% calculation of default edge orientation for each edge, in the
% direction of node1->node2 as given by enodesForEdge = [node1,node2]
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