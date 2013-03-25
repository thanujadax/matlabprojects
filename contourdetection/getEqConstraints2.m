function [Aeq,beq] = getEqConstraints2(numEdges,jEdges,edges2pixels)

[~, numJtypes] = size(jEdges);
% type 1 is J2 - junction with just 2 edges
nodeTypeStats = zeros(numJtypes,2);
% each row corresponds to a junction type. row 1: type 1 (J2)
% column 1: n.o. junction nodes of this type
% column 2: n.o. edge pair combinations to be activated
totJunctionVar = zeros(numJtypes,1); % stores the number of coefficients for each type of J
for i=1:numJtypes
    jEdges_i = jEdges{i};
    if(jEdges_i==0)
        % no junctions of this type
        nodeTypeStats(i,1) = 0;
        nodeTypeStats(i,2) = 0;
        totJunctionVar(i) = 0;
    else
        [numJ_i,numEdges_i] = size(jEdges_i);
        nodeTypeStats(i,1) = numJ_i;
        nodeTypeStats(i,2) = numEdges_i;
        numEdgeCombinations = nchoosek(numEdges_i,2);
        totJunctionVar(i) = numJ_i.*(numEdgeCombinations + 1); % 1 for the inactive node
        % clear nodeAngleCost_i
    end
end


% num cols of Aeq = 2*numEdges + sum_j(numNodes_j*numCombinations+1)
numCols_Aeq = 2*numEdges + sum(totJunctionVar); 
% num cols of Aeq = numEdges + 2*numJ
numRows_Aeq = numEdges + sum(nodeTypeStats(:,1))*2;

Aeq = zeros(numRows_Aeq,numCols_Aeq);
beq = zeros(numRows_Aeq,1);
beq(1:(numEdges + sum(nodeTypeStats(:,1)))) = 1;

%% activation/inactivation constraints for each edge
j = 1;
for i=1:numEdges
    Aeq(i,j:(j+1)) = 1;
    j = j+2;
end

%% activation/inactivation constraints for each junction node
colStop = numEdges*2;
rowStop = numEdges;
for jType = 1:numJtypes
    % for each junction type
    numNodes_j = nodeTypeStats(jType,1);
    numEdgePJ = jType + 1;      % number of edges per junction
    numCoef = nchoosek(numEdgePJ,2) + 1; % num edge pair combinations + inactivation  
    rowStart = rowStop + 1; 
    rowStop = rowStart - 1 + numNodes_j;
    for row=rowStart:rowStop
        colStart = colStop + 1;
        colStop = colStart - 1 + numCoef;
        Aeq(row,colStart:colStop) = 1;
    end    
end

% % J3
% j = numEdges*2+1;
% for i=(numEdges+1):(numEdges+numJ3)
%     Aeq(i,j:(j+3)) = 1;
%     j = j+4;
% end
% % J4
% j = numEdges*2+numJ3*4+1;
% for i=(numEdges+numJ3+1):(numEdges+numJ3+numJ4)
%     Aeq(i,j:(j+6)) = 1;
%     j = j+7;
% end

%% closedness constraints
% for each node, either exactly two edges should be active (active node) or
% all the edges should be inactive (inactive node)
numNodesTot = sum(nodeTypeStats(:,1));
rowStop = numEdges + numNodesTot;
jColIdStop = numEdges*2;
for jType=1:numJtypes
    % for each junction type
    % for each node, get the indices of the activeState edge variables
    jEdges_j = jEdges{jType};
    jEdgesOrderedInd_j = jEdges_j;
    % jEdges_j: the edge labels here follow the initial indexing. However,
    % the order of edges considered is slightly different to that due to
    % the removal of self edges and dangling edges. the edgeID->colInd of
    % the edge should be determined using the entries of the 1st column of
    % edges2pixels. This row number gives the edgeId to be used here.
    
    numNodes_j = nodeTypeStats(jType,1);
    if(numNodes_j~=0)
        for i=1:numNodes_j
            % for each node, get the edges fron jEdges_j
            % for each edge, get the row number from edges2pixels
            % add to jEdgesInd_j
            edgeIDs = jEdges_j(i,:);
            for m=1:numel(edgeIDs)
               edgeOrderInd = find(edges2pixels(:,1)==edgeIDs(m)); 
               jEdgesOrderedInd_j(i,m) = edgeOrderInd; 
            end


        end
    
    
        activeEdgeColInds = jEdgesOrderedInd_j.*2;     % edges are represented in pairs of state variables
                                            % the second element corresponds to
                                            % the active state
        % for each junction for type j
        rowStart = rowStop + 1;
        rowStop = rowStart - 1 + numNodes_j;
        k = 1;
        for row=rowStart:rowStop
            Aeq(row,activeEdgeColInds(k,:)) = 1; 
            k = k + 1;
            numEdges_i = nodeTypeStats(jType,2);
            numEdgeCombinations_j = nchoosek(numEdges_i,2);
            jColIdStart = jColIdStop + 2;
            jColIdStop = jColIdStart + numEdgeCombinations_j - 1;
            jIds = jColIdStart:jColIdStop;
            Aeq(row,jIds) = -2;               % refer closedness constraint formulation
            % jColIdStop = jColIdStop + numEdgeCombinations_j;
        end
    end
        
end



% % for all nodes, get the active edge state variable indices
% j3ActiveEdgeColInds = j3Edges .*2;
% j4ActiveEdgeColInds = j4Edges .*2;
% 
% % J3
% jColId = numEdges*2;
% k = 1;
% for i=(numEdges+numJ3+numJ4+1):(numEdges+numJ3*2+numJ4)
%    Aeq(i,j3ActiveEdgeColInds(k,:)) = 1;          % marks active edges for junction i
%    k = k+1;
%    jIds = (jColId+2):(jColId+4);        % activate the 3 active states coeff for J3
%    Aeq(i,jIds) = -2;                      % refer closedness constraint formulation
%    
%    % finally
%    jColId = jColId + 4;
% end
% 
% % J4
% jColId = numEdges*2 + numJ3*4;
% k = 1;
% for i=(numEdges+numJ3*2+numJ4+1):(numEdges+numJ3*2+numJ4*2)
%    Aeq(i,j4ActiveEdgeColInds(k,:)) = 1;          % marks active edges for junction i
%    k = k+1;
%    jIds = (jColId+2):(jColId+7);        % activate the 6 active states coeff for J4
%    Aeq(i,jIds) = -2;                      % refer closedness constraint formulation
%    
%    % finally
%    jColId = jColId + 7;
% end

