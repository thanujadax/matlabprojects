function f = getILPObjectiveVectorParametric(edgeUnary,nodeAngleCosts,...
            regionUnary,w_off_e,w_on_e,w_off_n,w_on_n,w_off_r,w_on_r)
numEdges = size(edgeUnary,1);

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

numRegions = numel(regionUnary);
numElements = numEdges*2 + sum(totJunctionVar) + numRegions*2;
f = zeros(numElements,1);
% order of elements in f
%[{edgeInactive},{edgeActive},{J2inactive},{J2Active},{J3inactive},{J4Active_3}...]

%% edge variables
j=1;
for i=1:2:2*numEdges
    f(i) = edgeUnary(j) * w_off_e;        % inactivation cost *w_off
    f(i+1) = edgeUnary(j) * w_on_e;     % activation cost for the same edge *w_on
    j = j+1;
end

f_stop_ind = numEdges*2;

%%  junction variables
for i=1:numJtypes
    % for each junction type
    clear nodeAngleCost_i
    nodeAngleCost_i = nodeAngleCosts{i};
    nodeAngleCost_i = nodeAngleCost_i .* w_on_n;
    [numJ,~] = size(nodeAngleCost_i);
    if(~isnan(nodeAngleCost_i))
        % inactivation cost - comes from learned param
        offJcost = ones(numJ,1).*w_off_n;
        
    % REMOVED: bbJunction cost assignment.
        
        nodeAngleCost_i = [offJcost nodeAngleCost_i];
        numCoeff_i = totJunctionVar(i);
        f_start_ind = f_stop_ind + 1;
        f_stop_ind = f_start_ind + numCoeff_i - 1;
        % assign coefficients to vector f
        angleCostMat_i = nodeAngleCost_i';
        f(f_start_ind:f_stop_ind) = angleCostMat_i(1:numCoeff_i);
    end
end

%% cells
k=1;
for i=(f_stop_ind+1):2:(f_stop_ind+numRegions*2)
    f(i) = regionUnary(k) * w_off_r; % inactivation
    f(i+1) = regionUnary(k) * w_on_r; % activation
    k = k + 1;
end
