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
    % listInds = junctionTypeListInds((junctionTypeListInds(:,jType)>0),jType);
    [rL,cL] = find(junctionTypeListInds(:,jType)>0);
    if(~isempty(rL))
        listInds = junctionTypeListInds(rL,jType);
        numColNodeEdges_i = jType+2;
        edgeSet = nodeEdges(listInds,2:numColNodeEdges_i);
        [r,c] = find(edgeSet>0);
        rmax = max(r);
        cmax = max(c);
        % edgeSetNoZeros = zeros(rmax,cmax);
        % edgeSetNoZeros(r,c) = edgeSet(r,c);
        % jEdges{jType} = edgeSetNoZeros;
        jEdges{jType} = edgeSet;
    else
        jEdges{jType} = 0;
    end
end