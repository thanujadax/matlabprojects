function outwardScores = getPreferredDirections(edges2nodes,junctionTypeListInds,...
            jAnglesAll,jAnglesAll_alpha,jEdges,edgeListInds)

% Inputs:
%   jEdges: cell array. each cell contains edgeIDs at each node for a
%   specific junctionType.

% Output: 
% outwardness: score reflecting the outwardness N1->N2 for the directional
% edges given in 
        
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
        [~,eLIDs] = edgeIDs_i(j,:);
        % eLID_directional = eLID or eLID+numEdges
        % get thetas
        thetas = thetas_i(j,:);
        % get alphas
        alphas = alphas_i(j,:);
        % get outwardness
        outwardness_j = getOutwardness(thetas,alphas);
        
        % get start node id
        N1_LID = junctionTypeListInds(j,i);
        % get end node ids for each edge
        nodes_eLIDs = edges2nodes(eLIDs,:);
        
        
        % get edgeLID_directional
        
    end
end