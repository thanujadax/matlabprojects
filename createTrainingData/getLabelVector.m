function labelVector = getLabelVector...
    (activeEdgeListInds,activeNodeListInds,activeRegionListInds,...
    numEdges,numRegions,jEdges,junctionTypeListInds,edgeListInds)

% create label vector
[~, numJtypes] = size(jEdges);
% type 1 is J2 - junction with just 2 edges
nodeTypeStats = zeros(numJtypes,2);
% each row corresponds to a junction type. row 1: type 1 (J2)
% column 1: n.o. junction nodes of this type
% column 2: n.o. number of node configurations (active and inactive)
totJunctionVar = zeros(numJtypes,1); % stores the number of coefficients for each type of J
totActiveJunctionConfs = zeros(numJtypes,1);

for i=1:numJtypes
    jEdges_i = jEdges{i};
    if(jEdges_i==0)
        % no junctions of this type
        nodeTypeStats(i,1) = 0;
        nodeTypeStats(i,2) = 0;
        totJunctionVar(i) = 0;
        totActiveJunctionConfs(i) = 0;
    else
        [numJ_i,numEdges_i] = size(jEdges_i);
        nodeTypeStats(i,1) = numJ_i;
        
        numEdgeCombinations = nchoosek(numEdges_i,2);
        nodeTypeStats(i,2) = numEdgeCombinations + 1; % +1 for the inactivation state
        
        totJunctionVar(i) = numJ_i.*(numEdgeCombinations + 1); % 1 for the inactive state
        totActiveJunctionConfs(i) = numJ_i.*numEdgeCombinations;
        % clear nodeAngleCost_i
    end
end

numElements = numEdges*2 + sum(totJunctionVar) + numRegions*2;

labelVector = zeros(numElements,1);

%% Edge activation
edgeListInd=1;
for i=1:2:2*numEdges
    if(sum(ismember(activeEdgeListInds,edgeListInd))>0)
        % edge is active
        labelVector(i+1) = 1;
    else
        % edge is inactive
        labelVector(i) = 1;
    end
    
    edgeListInd = edgeListInd + 1;
end

%% Node activation
% given the active edge list inds, which node configurations are active?
f_stop_ind = 2*numEdges;
nodeInd_x = 0;
for i=1:numJtypes
    % for each junction type
    clear nodeAngleCost_i
    nodeEdgeIDs_i = jEdges{i};
    numNodes_i = size(nodeEdgeIDs_i,1);
    numCoeff_i = nodeTypeStats(i,2);
    if( sum(sum(~isnan(nodeEdgeIDs_i)))>0 && nodeEdgeIDs_i(1)~=0)
        for j=1:numNodes_i
            % check if the node is an active node as labeled
            nodeCoeff = zeros(1,numCoeff_i);
            nodeInd_x = nodeInd_x + 1;
            nodeListInd = junctionTypeListInds(j,i);
            if(sum(ismember(activeNodeListInds,nodeListInd))>0)
                % node is labeled active
                nodeLabel = 1;  % flag
                [~, nodeEdgeListInds_i] = intersect(edgeListInds,nodeEdgeIDs_i(j,:));
                activeEdgeLI_logical = ismember(activeEdgeListInds,nodeEdgeListInds_i);
                activeEdgeListInds_j = activeEdgeListInds(activeEdgeLI_logical);
                numActiveEdges_i = sum(activeEdgeLI_logical);
                % if all the edges connected to the node are inactive, the node is
                % inactive
                if(numActiveEdges_i<2)
                    nodeCoeff(1) = 1;
                elseif(numActiveEdges_i==2)
                    nodeCoeff(1) = 0;
                    % find which active configuration
                    activeConfInd = getNodeActiveConfig...
                        (nodeEdgeListInds_i,activeEdgeListInds_j);
                    if(activeConfInd>0)
                        nodeCoeff(activeConfInd+1) = 1;
                    else
                        disp('**********************************************')
                        disp('ERROR1: getLabelVector')
                    end
                elseif(numActiveEdges_i>2)
                    disp('**********************************************')
                    disp('ERROR2: getLabelVector. more than 2 edges active!')
                    
                else
                    disp('**********************************************')
                    disp('ERROR3: getLabelVector. active node config unresolvable')
                end


            else
                % node is labeled inactive
                nodeLabel = 0; % flag
                % mark inactive
                nodeCoeff(1) = 1;
            end

            f_start_ind = f_stop_ind + 1;
            f_stop_ind = f_start_ind + numCoeff_i - 1;
            % assign coefficients to vector f

            labelVector(f_start_ind:f_stop_ind) = nodeCoeff(1:numCoeff_i);

        end
    end
end
%% Region activation
regionListInd=1;
f_start_ind = f_stop_ind + 1;
f_stop_ind = f_stop_ind + numRegions*2;
for i=f_start_ind:2:f_stop_ind
    if(sum(ismember(activeRegionListInds,regionListInd))>0)
        % edge is active
        labelVector(i+1) = 1;
    else
        % edge is inactive
        labelVector(i) = 1;
    end
    
    regionListInd = regionListInd + 1;
end

%% supplementary functions

function activeConfInd = getNodeActiveConfig...
                        (nodeEdges_i,activeEdgeListInds_i)
activeConfInd = 0;
numEdgesPerNode = numel(nodeEdges_i);
edgeLIDvect = 1:numEdgesPerNode;
combinations = nchoosek(edgeLIDvect,2);
numCombinations = size(combinations,1);
for i=1:numCombinations
   combinationLID_i = combinations(i,:);
   combination_i = nodeEdges_i(combinationLID_i);
   if(sum(ismember(combination_i,activeEdgeListInds_i))==2)
       % right combination!
       activeConfInd = i;
       break
   end
end

