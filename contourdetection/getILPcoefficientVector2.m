function f = getILPcoefficientVector2(edgePriors,nodeAngleCosts)
numEdges = size(edgePriors,1);

[~, numJtypes] = size(nodeAngleCosts);
% type 1 is J2 - junction with just 2 edges
nodeTypeStats = zeros(numJtypes,2);
% each row corresponds to a junction type. row 1: type 1 (J2)
% column 1: n.o. junction nodes of this type
% column 2: n.o. edge pair combinations to be activated
totJunctionVar = zeros(numJtypes,1); % stores the number of coefficients for each type of J
for i=1:numJtypes
    nodeAngleCost_i = cell2mat(nodeAngleCosts{i});
    [numJ_i,numCombinations] = size(nodeAngleCost_i);
    nodeTypeStats(i,1) = numJ_i;
    nodeTypeStats(i,2) = numCombinations;
    totJunctionVar(i) = nodeTypeStats(i,1).*nodeTypeStats(i,2);
    clear nodeAngleCost_i
end


numElements = numEdges*2 + sum(totJunctionVar);
f = zeros(numElements,1);
% order of elements in f
%[{edgeInactive},{edgeActive},{J2inactive},{J2Active},{J3inactive},{J4Active_3}...]

% edge variables
j=1;
for i=1:2:2*numEdges
    f(i) = edgePriors(j);        % inactivation cost
    f(i+1) = -edgePriors(j);     % activation cost for the same edge
    j = j+1;
end

f_stop_ind = numEdges*2;
% junction variables
for i=1:numJtypes
    % for each junction type
    nodeAngleCost_i = cell2mat(nodeAngleCosts{i});
    maxJcost = max(nodeAngleCost_i,[],2);          % inactivation cost
    nodeAngleCost_i = [maxJcost nodeAngleCost_i];
    numCoeff_i = totJunctionVar(i);
    f_start_ind = f_stop_ind + 1;
    f_stop_ind = f_start_ind + numCoeff_i - 1;
    % assign coefficients to vector f
    f(f_start_ind:f_stop_ind) = nodeAngleCost_i(1:numCoeff_i);
end

