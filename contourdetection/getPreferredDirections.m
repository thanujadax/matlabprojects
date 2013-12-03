function outwardness = getPreferredDirections(edges2nodes,junctionTypeListInds,...
            jAnglesAll,jAnglesAll_alpha,jEdges,edgeListInds)

% Inputs:
%   jEdges: cell array. each cell contains edgeIDs at each node for a
%   specific junctionType.
        
numEdges = size(edges2nodes,1);

outwardness = zeros(numEdges,1);

[~,numJtypes] = size(junctionTypeListInds);

for i=1:numJtypes
    numJ_i = sum(junctionTypeListInds(:,i)>0);
    edgeIDs_i = jEdges{i};
    thetas_i = jAnglesAll{i};
    alphas_i = jAnglesAll_alpha{i};
    
    for j=1:numJ_i
        
    end
end