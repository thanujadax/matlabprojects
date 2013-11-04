function qsparse = getQuadraticObjective_PE(edgePriors,nodeAngleCosts,...
            regionPriors,numParam)
        
numEdges = size(edgePriors,1);

[~, numJtypes] = size(nodeAngleCosts);
% type 1 is J2 - junction with just 2 edges
nodeTypeStats = zeros(numJtypes,2);
% each row corresponds to a junction type. row 1: type 1 (J2)
% column 1: n.o. junction nodes of this type
% column 2: n.o. edge pair combinations to be activated
totJunctionVar = zeros(numJtypes,1); % stores the number of coefficients for each type of J
for i=1:numJtypes
    nodeAngleCost_i = nodeAngleCosts{i};
    if(isnan(nodeAngleCost_i))
        % ignore - no such junctions of this type
        nodeTypeStats(i,1) = 0;
        nodeTypeStats(i,2) = 0; 
        totJunctionVar(i) = 0;
    else
        [numJ_i,numCombinations] = size(nodeAngleCost_i);
        numCombinations = numCombinations + 1;  % 1 for the inactive junction
        nodeTypeStats(i,1) = numJ_i;
        nodeTypeStats(i,2) = numCombinations; 
        totJunctionVar(i) = nodeTypeStats(i,1).*nodeTypeStats(i,2);
        clear nodeAngleCost_i
    end
end

% numParam = 7;
numRegions = numel(regionPriors);
numIntVariables = numEdges*2 + sum(totJunctionVar) + numRegions*2;
totNumVariables = numIntVariables + numParam;
qmat = zeros(numIntVariables,numParam);
% order of elements in Q
%[{edgeInactive},{edgeActive},{J2inactive},{J2Active},{J3inactive},{J4Active_3}...,{RegionInact},{RegionActive}]

colOffset = 0;

%% edge variables
j=1;
for i=1:2:2*numEdges
    % col1: inactivation weight
    colID = colOffset + 1;
    qmat(i,colID) = edgePriors(j);        % inactivation cost *w_off
    % col2: activation weight
    colID = colOffset + 2;
    qmat(i+1,colID) = -edgePriors(j);     % activation cost for the same edge *w_on
    j = j+1;
end

q_stop_ind = numEdges*2;

%%  junction variables

for i=1:numJtypes
    % for each junction type
    clear nodeAngleCost_i
    nodeAngleCost_i = nodeAngleCosts{i};
    if(~isnan(nodeAngleCost_i))
        % col 3: inactivation
        offJcost = 1; % correspond to a constant weight parameter to be inferred by QP
        q_stop_ind = q_stop_ind + 1;
        colID = colOffset + 3;
        qmat(q_stop_ind,colID) = offJcost;
        
        numJunctVar = numel(nodeAngleCost_i);
        for m=1:numJunctVar
            q_stop_ind = q_stop_ind + 1;
            junctCost = nodeAngleCost_i(m);
            if(junctCost>0)
                % positive jn cost: col 4
                colID = colOffset + 4;
                qmat(q_stop_ind,colID) = junctCost;
            else
                % negative or zero jn cost: col 5
                colID = colOffset + 5;
                qmat(q_stop_ind,colID) = junctCost;
            end
        end
    end
end

%% cells
k=1;
for i=(q_stop_ind+1):2:(q_stop_ind+numRegions*2)
    % col 6: region inactivation
    colID = colOffset + 6;
    qmat(i,colID) = regionPriors(k); % inactivation
    % col 7: region activation
    colID = colOffset + 7;
    qmat(i+1,colID) = - regionPriors(k); % activation
    k = k + 1;
end

%% Creating sparse output matrix
[r,c] = find(qmat);
c = c + numIntVariables;
s = qmat(qmat~=0);

qsparse = sparse(r,c,s,totNumVariables,totNumVariables);