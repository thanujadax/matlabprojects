function [inEdgeListInds, inEdgeIDs] = getInEdges(twoCellEdges,cellActivationVector,...
                edgeActivationVector,edges2cells,edgeIDsAll)

% Output:
%   inEdges: the vector containing edgeListIDs that are inactive, but located
%   inbetween two active regions

% Inputs:
%   twoCellEdges - edgeIDs which are in between 2 cells
%   cellActivationVector - contains 1 for cellIDs that are active, 0
%   otherwise
%   edgeActivationVector - contains 1 for edgeListInds  that are active, 0
%   otherwise
%   edges2cells - 
%   edgeIDsAll - list of all edgeIDs available in order

inEdgeListInds = [];
inEdgeIDs = [];

for i=1:numel(twoCellEdges)
    edgeListInd_i = find(edgeIDsAll==twoCellEdges(i));
    if(~edgeActivationVector(edgeListInd_i))
        % edge is inactive
        cellInd_1 = edges2cells(twoCellEdges(i),1);
        cellInd_2 = edges2cells(twoCellEdges(i),2);
        
        cellState_1 = cellActivationVector(cellInd_1);
        cellState_2 = cellActivationVector(cellInd_2);
        if(cellState_1>0 && cellState_2>0)
            % append this edgeListInd to the set of inEdges
            inEdgeListInds = [inEdgeListInds; edgeListInd_i];
            inEdgeIDs = [inEdgeIDs; twoCellEdges(i)];
        end
    end
end