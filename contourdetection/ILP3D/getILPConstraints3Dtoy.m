function [A,b,senseArray] = getILPConstraints3Dtoy...
            (c_edgeIDs,c_Ts,sectVarStats,c_TsetsWithCommonEdges)

% 2014.02.17

% Outputs:
%   A - sparse matrix containing coefficients for ILP
%   [sect1:{edges}{nodes}{regions}][sect2:...][]...[Ts_1-2][Ts_2-3]...
%   b - 
%   senseArray - 

% Inputs:
%   c_Ts - cell array containing Ts for each adjacent section pair. 
%       T_i = [sect1][edgeID][sect2][nodeID1][nodeID2]
%   Data struct info
%   Each T_i = c_Ts{i} has the column structure that should be used in the
%   code 
%   sectVarStats - contains the number of edge, node, region vars for each
%       section. rowID=sectID. cols: {edgeVars}{nodeConfs}{regionVars}
%   c_TsetsWithCommonEdges: cell array. Each cell contains a vector of
%   T_ids (including inactive state which is usually T_id + tot_Ts)

sect1ID_col = 1;
sect1EdgeID_col = 2;
sect2ID_col = 3;
sect2NodeID1_col = 4;
sect2NodeID2_col = 5;

%% Init
rowStop = 0;
numSections = size(c_edgeIDs,2);
numSectionPairs = size(c_Ts,2);
offSetTInds = sum(sum(sectVarStats)); % tot num of section vars {e,n,r}

A = sparse(totNumConstraints,totColsA);
b = zeros(totNumConstraints,1);
senseArray(1:totNumConstraints) = '=';
%% 2D constraints
% formulated for each of the sections

%% 3D constraints
% formulated for each adjacent pair of sections
%% T activation - binary constraint
% 
%% T activation - no overlapping Ts
% for each pair of sections, get the set of Ts based on each edge.

for secPair = 1:numSectionPairs
    Ts_i = c_Ts{secPair};
    numEdgesInSecPair = size(Ts_i,1);
    
    for edge_LInd = 1:numEdgesInSecPair
        sectIDforEdge_i = Ts_i(edge_LInd,sect1ID_col);
        edgeID_in_sect_i = Ts_i(edge_LInd,sect1EdgeID_col); 
        
        tInds_edge_i = getTIndsForEdgeGivenSectionPair...
                (Ts_i,edgeID_in_sect_i,sectIDforEdge_i); 
        % does not include t_offStates
        
        tInds_edge_i_offset = tInds_edge_i + offSetTInds;
        rowStop = rowStop + 1;
        A(rowStop,tInds_edge_i_offset) = 1.1;
        senseArray(rowStop) = '<'; 
        b(rowStop) = 1;
    end
end

%% T activation - activating relevant edges and nodes on the adjacent slices


%% T - Laterally contiguous
% For each active T, there should be another active T on either side of it,
% sticking to its sides in the z direction

% get sets of Ts where out of each set only 2 or 0 can be active. The first
% element in each set stands for the inactive state
numTsets = size(c_TsetsWithCommonEdges,2);
for i=1:numTsets
    Tset_i = c_TsetsWithCommonEdges{i};
    inactiveT_id = Tset_i(1);
    Tset_i(1) = []; % now contains active T_ids only
    
    colID_inactive_state = inactiveT_id + offSetTInds;
    colIDs_active_states = Tset_i + offSetTInds;
    
    rowStop = rowStop + 1;
    A(rowStop,colID_inactive_state) = 2;
    A(rowStop,colIDs_active_states) = 1;
    b(rowStop) = 2;
    % senseArray(rowStop) = '='; Default value 
end

