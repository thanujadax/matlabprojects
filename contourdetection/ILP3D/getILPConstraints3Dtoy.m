function [A,b,senseArray] = getILPConstraints3Dtoy(c_edgeIDs,c_Ts)

% 2014.02.17

% Outputs:
%   A - sparse matrix containing coefficients for ILP
%   b - 
%   senseArray -

% Inputs:
%   c_Ts - cell array containing Ts for each adjacent section pair. 
%       T_i = [sect1][edgeID][sect2][nodeID1][nodeID2]

% Data struct info
% Each T_i = c_Ts{i} has the column structure that should be used in the
% code 
sect1ID_col = 1;
sect1EdgeID_col = 2;
sect2ID_col = 3;
sect2NodeID1_col = 4;
sect2NodeID2_col = 5;

%% Init
rowStop = 0;
numSections = size(c_edgeIDs,2);
numSectionPairs = size(c_Ts,2);

A = sparse(totRowsA,totColsA);
%% 2D constraints
% formulated for each of the sections

%% 3D constraints
% formulated for each adjacent pair of sections
%% T activation - no overlapping Ts
% for each pair of sections, get the set of Ts based on each edge.

for secPair = 1:numSectionPairs
    Ts_i = c_Ts{secPair};
    numEdgesInSecPair = size(Ts_i,1);
    
    for edge_LInd = 1:numEdgesInSecPair
        sectIDforEdge_i = Ts_i(edge_LInd,sect1ID_col);
        edgeID_in_sect_i = Ts_i(edge_LInd,sect1EdgeID_col); 
        tInds_edge_i = getTIndsForEdgeGivenSectionPair(); % includes t_offState
        tInds_edge_i_offset = tInds_edge_i + offSetTInds;
        rowStop = rowStop + 1;
        A(rowStop,tInds_edge_i_offset) = 1;
        senseArray(rowStop) = '=';
        b(rowStop) = 1;
    end
end
%% T - Laterally contiguous

%% T activation - activating relevant edges and nodes on the adjacent slices

