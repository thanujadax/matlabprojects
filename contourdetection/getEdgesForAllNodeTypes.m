function jEdges = getEdgesForAllNodeTypes(nodeEdges,junctionTypeListInds)
% returns a cell array.
% jEdges{i} - each row corresponds to the set of edges for each
% junction of type i (type1 = J2)

% Inputs:
%   nodeEdges -
%   junctionTypeListInds - each column contains the list indices of
%   the nodes of a junction type wrt the order given in nodeEdges.

[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);
maxNumEdgesPerNode = numJtypes + 1;     % list starts with J2 (= two edges)
[numNodes,numColNodeEdges] = size(nodeEdges);

% jEdges = zeros(maxNodesPerJtype,maxNumEdgesPerNode,numJtypes);
jEdges = cell(1,numJtypes);

for jType=1:numJtypes
    listInds = junctionTypeListInds((junctionTypeListInds(:,jType)>0),jType);
    if(~isempty(listInds))
        edgeSet = nodeEdges(listInds,2:numColNodeEdges);
        [r c] = find(edgeSet>0);
        edgeSetNoZeros(r,c) = edgeSet(r,c);
        jEdges{jType} = edgeSetNoZeros;
    else
        jEdges{jType} = 0;
    end
end