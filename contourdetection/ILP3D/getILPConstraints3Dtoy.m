function [A,b,senseArray,numEdges,numNodeConf,numRegions,nodeTypeStats] = ...
            getILPConstraints3Dtoy...
            (c_edgeIDs,c_Ts_on,sectVarStats,c_TsetsWithCommonEdges,...
            edgeListInds,edges2nodes,nodeEdgeIDs,junctionTypeListInds,...
            jEdges,dirEdges2regionsOnOff,setOfRegions,activeWSregionListInds_tr)

% 2014.02.17

% Outputs:
%   A - sparse matrix containing coefficients for ILP
%   [sect1:{edges}{nodes}{regions}][sect2:...][]...
%       [TsOn_1-2][TsOn_2-3]...[TsOff_1-2][TsOff_2-3]
%   b - 
%   senseArray - 

% Inputs:
%   c_Ts_on - cell array containing Ts for each adjacent section pair. 
%       T_i = [sect1][edgeID][sect2][nodeID1][nodeID2]
%       Considers only the on state of the Ts. 
%   sectVarStats - contains the number of edge, node, region vars for each
%       section. rowID=sectID. cols: {edgeVars}{nodeConfs}{regionVars}
%   c_TsetsWithCommonEdges: cell array containing 3 cells.
%       First cell contains a scaler: T_id
%       Each of the other 2 cells in the cell array contains a vector of T_ids which 
%       stick to one side of the T in concern.


%   Data struct info
%   Each T_i = c_Ts_on{i} has the column structure that should be used in the
%   code
sect1ID_col = 1;
sect1EdgeID_col = 2;
sect2ID_col = 3;
sect2NodeID_col = 4;
sect1NodeID1_col = 5;
sect1NodeID2_col = 6;

%% Enable/disable constraints
% 2D

% 3D
activation_T = 1;
nonoverlap_T = 1;
T_nodes_edges = 1;
T_lateral_contig = 0;

%% 2D constraints
% formulated for each of the sections
[A2d,b2d,senseArray2d,numEdges,numNodeConf,numRegions,nodeTypeStats]...
    = getILPConstraints(edgeListInds,edges2nodes,nodeEdgeIDs,junctionTypeListInds,...
        jEdges,dirEdges2regionsOnOff,setOfRegions,activeWSregionListInds_tr);

%% Init
rowStop = 0;
numSections = size(c_edgeIDs,2);
numSectionPairs = size(c_Ts_on,2);

% stats for Ts
% table containing number of Ts_on for each adjacent section pair
numTsOn_sectPairs = zeros(numSectionPairs,1);
for secPair = 1:numSectionPairs
    Ts_on_i = c_Ts_on{secPair};
    numTsOn_sectPairs(secPair) = size(Ts_on_i,1);
end
numTs_on = sum(numTsOn_sectPairs);
offSetTInds_on = sum(sum(sectVarStats)); % tot num of section vars {e,n,r}
offSetTInds_off = offSetTInds_on + numTs_on;

A3d = sparse(totNumConstraints,totColsA);
b3d = zeros(totNumConstraints,1);
senseArray3d(1:totNumConstraints) = '=';


%% 3D constraints
% formulated for each adjacent pair of sections
%% T activation - binary constraint
if(activation_T)
    % tOn_i + tOff_i = 1;
    tmpNumTsDone = 0;
    for secPair = 1:numSectionPairs
        numTsOn_pair_i = size(c_Ts_on{secPair},1);
        % get the indices of each pair of on_off Ts in the variable vector.
        % enforce binary activation constraint
        Ts_on_seq_i = 1:numTsOn_pair_i;
        Ts_off_listInds = Ts_on_seq_i + offSetTInds_off + tmpNumTsDone;
        Ts_on_listInds = Ts_on_seq_i + offSetTInds_on + tmpNumTsDone;
        tmpNumTsDone = tmpNumTsDone + numTsOn_pair_i;

        % set values in A, senseArray and b.
        rowStart = rowStop + 1;
        rowStop = rowStop + numTsOn_pair_i;
        for i = rowStart:rowStop
            col_activeState = Ts_on_listInds(i);
            col_inactiveState = Ts_off_listInds(i);
            A3d(i,col_activeState) = 1;
            A3d(i,col_inactiveState) = 1;
            b3d(i) = 1;
            % senseArray(i) = '='; % default value
        end

    end
end
%% T activation - no overlapping Ts
if(nonoverlap_T)
    % for each pair of sections, get the set of Ts based on each edge.
    tmpNumTsDone = 0;
    for secPair = 1:numSectionPairs
        Ts_on_i = c_Ts_on{secPair};
        numEdgesInSecPair = size(Ts_on_i,1);

        for edge_LInd = 1:numEdgesInSecPair
            sectIDforEdge_i = Ts_on_i(edge_LInd,sect1ID_col);
            edgeID_in_sect_i = Ts_on_i(edge_LInd,sect1EdgeID_col); 

            tInds_edge_i = getTIndsForEdgeGivenSectionPair...
                    (Ts_on_i,edgeID_in_sect_i,sectIDforEdge_i); 
            % does not include t_offStates

            tInds_edge_i_offset = tInds_edge_i + offSetTInds_on + tmpNumTsDone;
            rowStop = rowStop + 1;
            A3d(rowStop,tInds_edge_i_offset) = 1.1;
            senseArray3d(rowStop) = '<'; 
            b3d(rowStop) = 1;
        end
        tmpNumTsDone = tmpNumTsDone + size(Ts_on_i,1);
    end
end

%% T activation - activating relevant edges and nodes on the adjacent slices
if(T_nodes_edges)
    % e + n + T_on -3T_off <=0;

    tmpNumTsDone = 0;
    % for each T
    for secPair = 1:numSectionPairs
        Ts_on_i = c_Ts_on{secPair};
        numTs_i = size(Ts_on_i);

        % get the list of section ids
        sect1IDs_i = Ts_on_i(:,sect1ID_col); % contains edge of T
%         sect2IDs_i = Ts_on_i(:,sect2ID_col); % contains node
%         % get set of nodes
%         nodeLIDs_i = Ts_on_i(:,sect2NodeID_col);
%         nodeColIDs_i = get
        % get set of T_on_ids
        T_on_colIDs = (1:numTs_i) + offSetTInds_on + tmpNumTsDone;
        % get set of T_off_ids
        T_off_colIDs = (1:numTs_i) + offSetTInds_off + tmpNumTsDone;

        % get set of edges
        edgeIDs_i = Ts_on_i(:,sect1EdgeID_col);

        [edgeID_cols_i,nodeID_cols_i] = getEdgeNodeIDCols...
                    (edgeIDs_i,nodeIDs_i,sect1IDs_i,sectVarStats);

        for i=1:numTs_i
            % for each T, enforce e + n + T_on - 3T_off <=0

            % set A,b,senseArray
            rowStop = rowStop + 1;
            colIDs_1 = [edgeID_cols_i(i) nodeID_cols_i(i) T_on_colIDs(i)];
            A3d(rowStop,colIDs_1) = 1;
            A3d(rowStop,T_off_colIDs(i)) = -3;
            b3d = 0.1;
            senseArray3d = '<';

        end

        tmpNumTsDone = tmpNumTsDone + numTs_i; 
    end
end

%% T - Laterally contiguous
if(T_lateral_contig)
    % For each active T, there should be another active T on either side of it,
    % sticking to its sides in the z direction

    % get sets of Ts where out of each set only 2 or 0 can be active. The first
    % element in each set stands for the inactive state
    numTsets = size(c_TsetsWithCommonEdges,2);
    for i=1:numTsets
        Tset_i = c_TsetsWithCommonEdges{i}; % 2 vectors
        inactiveT_id = Tset_i(1);
        Tset_i(1) = []; % now contains active T_ids only

        colID_inactive_state = inactiveT_id + offSetTInds_on;
        colIDs_active_states = Tset_i + offSetTInds_on;

        rowStop = rowStop + 1;
        A3d(rowStop,colID_inactive_state) = 2;
        A3d(rowStop,colIDs_active_states) = 1;
        b3d(rowStop) = 2;
        % senseArray(rowStop) = '='; Default value 
    end
end

%% Concatanate 2D and 3D matrices
% A
[is, js, ss] = find(A2d);
[it, jt, st] = find(A3d);
A = sparse([is; it + size(S,1)], [js; jt], [ss; st]);
b = [b2d; b3d];
senseArray = [senseArray2d senseArray3d];
