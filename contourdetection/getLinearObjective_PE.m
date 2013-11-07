function f = getLinearObjective_PE(edgePriors,nodeAngleCosts,...
            bbJunctionsListInds,junctionTypeListInds,bbJunctionCost,...
            regionPriors,activeEdgeListInds,activeRegionListInds,...
            activeNodeListInds,edgeReward,regionReward,nodeReward,numParam)
        
% linear objective function to be minimized in parameter estimation

% Inputs:
%   regionReward - negative value (decreases cost)
        
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

numRegions = numel(regionPriors);
numElements = numEdges*2 + sum(totJunctionVar) + numRegions*2;
f = zeros(numElements,1);
% order of elements in f
%[{edgeInactive},{edgeActive},{J2inactive},{J2Active},{J3inactive},{J4Active_3}...]

%% edge variables
j=1;
for i=1:2:2*numEdges
    if(sum(ismember(activeEdgeListInds,j))>0)
        % active edge
        f(i) = -edgeReward;     % penalize
        f(i+1) = edgeReward;    % reward
    else
        % inactive edge
        f(i) = edgeReward;     % reward
        f(i+1) = - edgeReward;    % penalize
    end
    j = j+1;
end

f_stop_ind = numEdges*2;

%%  junction variables
for i=1:numJtypes
    % for each junction type
    clear nodeAngleCost_i
    nodeAngleCost_i = nodeAngleCosts{i};untitled10.m
    if(~isnan(nodeAngleCost_i))
        maxJcost = max(nodeAngleCost_i,[],2);          % inactivation cost
        minJcost = min(nodeAngleCost_i,[],2);          
%         offJcost = (minJcost + maxJcost)/2;
        % offJcost = -1 .* minJcost;
%         offJcost = 0;
        offJcost = 1;  % will have to be weighted by a parameter that has to be learned.
        
        % identify the bbJunctions and set a very high inactivation cost
        clear nodeListInds_i
        nodeListInds_i = junctionTypeListInds(:,i); % list inds of nodes of type i
        [bbJunctionListInds_i,bbListInds_i,~] = intersect(nodeListInds_i,bbJunctionsListInds);
        if(numel(bbJunctionListInds_i)>0)
            % assign a very high inactivation cost for this nodes
            % (avgJcost)
            offJcost(bbListInds_i) = bbJunctionCost;
        end
        
    %     nodeAngleCost_i = [maxJcost nodeAngleCost_i];
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
rID=0;
for i=(f_stop_ind+1):2:(f_stop_ind+numRegions*2)
   rID = rID + 1;
   if(sum(ismember(activeRegionListInds,rID))>0)
       % label: active
       f(i) = -regionReward; % inactive state, cost +
       f(i+1) = regionReward; % active state, reward -
   else
       % label: inactive
       f(i) = regionReward; % inactive state
       f(i+1) = -regionReward; % active state
   end
    
end
fStop = i+1;
%%  params to optimize
fStart = fStop + 1;
fStop = fStop + numParam;
f(fStart:fStop) = 0;

