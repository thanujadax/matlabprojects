function [outwardScores,edges2nodes_directional] = getPreferredDirections...
                (edges2nodes,junctionTypeListInds,jAnglesAll,...
                jAnglesAll_alpha,jEdges,edgeListInds)

% Inputs:
%   jEdges: cell array. each cell contains edgeIDs at each node for a
%   specific junctionType.

% Output: 
% outwardness: score reflecting the outwardness N1->N2 for the directional
% edges given in edges2nodes_directional
        
numEdges = size(edges2nodes,1);

edges2nodes_complement = edges2nodes;
edges2nodes_complement(:,1) = edges2nodes(:,2);
edges2nodes_complement(:,2) = edges2nodes(:,1);
edges2nodes_directional = [edges2nodes; edges2nodes_complement];

outwardScores = zeros(numEdges*2,1);

[~,numJtypes] = size(junctionTypeListInds);

for i=1:numJtypes
    numJ_i = sum(junctionTypeListInds(:,i)>0);
    edgeIDs_i = jEdges{i};
    thetas_i = jAnglesAll{i};
    alphas_i = jAnglesAll_alpha{i};
    numEdgesPerNode = size(edgeIDs_i,2);
    
    for j=1:numJ_i
        % get edgeLIDs
        [~,eLIDs] = intersect(edgeListInds,edgeIDs_i(j,:));
        % eLID_directional = eLID or eLID+numEdges
        % get thetas
        thetas = thetas_i(j,:);
        % get alphas
        alphas = alphas_i(j,:);
        % get outwardness
        outwardness_j = getOutwardness(thetas,alphas);
        
        % get start node id (this node id)
        N1_LID = junctionTypeListInds(j,i);
        % get end node ids for each edge
        nodes_all_eLIDs = edges2nodes(eLIDs,:); % n-by-2 matrix
        % some of the edges will have this node as N1. for these edges,
        % eLIDs should not be adjusted. For the rest, 
        % edgeLIDs = edgeLIDs + numEdges
        
        
        % eLIDs directional
        eLID_N1N2_logicalPos_eLIDs = logical(nodes_all_eLIDs(:,1)==N1_LID);
        eLID_N2N1_logicalPos_eLIDs = logical(nodes_all_eLIDs(:,2)==N1_LID);
        eLID_N1N2 = eLIDs(eLID_N1N2_logicalPos_eLIDs);
        eLID_N2N1 = eLIDs(eLID_N2N1_logicalPos_eLIDs);
        eLID_N2N1 = eLID_N2N1 + numEdges; % shifted in directed edge index
        
        outwardScores(eLID_N1N2) = outwardness_j(eLID_N1N2_logicalPos_eLIDs);
        outwardScores(eLID_N2N1) = outwardness_j(eLID_N2N1_logicalPos_eLIDs);
        
    end
end