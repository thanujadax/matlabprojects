function [jEdges,jAngles,dTheta] = getInfoAllNodeTypes(nodeEdges,junctionTypeListInds,...
            nodeInds,edges2pixels,orientedScoreSpace3D,sizeR,sizeC,angleStep)
% returns a 3D array.
% jEdges(:,:,1) - each row corresponds to the set of edges for each
% junction of type 1 (= J2)

% Inputs:
%   nodeEdges -
%   junctionTypeListInds - each column contains the list indices of
%   the nodes of a junction type wrt the order given in nodeEdges.

[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);
maxNumEdgesPerNode = numJtypes + 1;     % list starts with J2 (= two edges)
[numNodes,numColNodeEdges] = size(nodeEdges);

jEdges = zeros(maxNodesPerJtype,maxNumEdgesPerNode,numJtypes);
jAngles = zeros(maxNodesPerJtype,maxNumEdgesPerNode,numJtypes);


for jType=1:numJtypes
    listInds = junctionTypeListInds(junctionTypeListInds(:,jType)>0);
    edgeSet = nodeEdges(listInds,2:numColNodeEdges);
    edgeSetNoZeros = edgeSet(edgeSet>0);
    [rows cols] = size(edgeSetNoZeros);
    % for each edge in the edgeSet of this type of node Jn
    for i=1:rows
        for j=1:cols
            jEdges(i,j,jType) = edgeSetNoZeros(i,j);
            % getAngles
            edgePixelInds = edges2pixels(edgeSetNoZeros(i,j),:);
            edgePixelInds = edgePixelInds(edgePixelInds>0);
            
        end
    end
    
end