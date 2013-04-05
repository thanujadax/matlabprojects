function junctionTypeListInds = getJunctionTypeListInds(nodeEdges)
% Input:
%   nodeEdges - first column contains the indices of all the junction
%   nodes. the rest of the entries in each row contains edgeIDs corresponding to 
%   each junction node.

% Output:
%   junctionTypeListInds -
%       row1- J2 junction node list index of nodeEdges
%       row2- J3 junction node list index of nodeEdges
%       ...

[numJunctions numMaxEdgesPJn] = size(nodeEdges);
numJunctionTypes = numMaxEdgesPJn - 2;

for i=1:numJunctions
    % get the number of edges at each junction i
    numEdgesAtJunction(i,1) = sum(nodeEdges(i,:)>0) -1;
end

for i=1:numJunctionTypes
    % for each junction type J2, J3 etc
    j = i+1;        % num of edges at junction type i. J2 (type 1) has 2.
    % get the nodeListIDs with j number of columns in nodeEdges
    junctionListInds_j = find(numEdgesAtJunction==j);
    if(~isempty(junctionListInds_j))
        for k=1:numel(junctionListInds_j)
            junctionTypeListInds(k,i) = junctionListInds_j(k);
        end
    else
        junctionTypeListInds(1,i) = 0; % in case a certain junction type has no instances 
    end
end


