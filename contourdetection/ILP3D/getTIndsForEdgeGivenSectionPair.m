function tIndsForEdge = getTIndsForEdgeGivenSectionPair(Ts,edgeID,sectID)

% Inputs:
%   Ts - matrix containing Ts for the ith adjacent section pair
%       columns of Ts_i are organized as follows
%       Ts_i = [sect1][edgeID][sect2][nodeID1][nodeID2]
%   edgeID - edgeID for which we should find all the Ts which are based on
%       this edge

% Output:
%   tIndsForEdge - T inds as given by the rowID of Ts_i

% get all row IDs for sect1(col1) = sectID
rowIDsForSectionID = find(Ts(:,1) == sectID);
% out of those rowIDs, pick the ones that have col2==edgeID.
TsForSection1 = Ts(rowIDsForSectionID,:);
rowIDsForEdgeIDs_logical = (TsForSection1(:,2)==edgeID);

tIndsForEdge = rowIDsForSectionID(rowIDsForEdgeIDs_logical);

