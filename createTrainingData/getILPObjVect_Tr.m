function f = getILPObjVect_Tr(labelImage,ws,edgeListInds,...
                setOfRegions,edges2nodes,numEdges,numNodeConf,numRegions,...
                edgeUnary)

% objective to be minimized (costs are positive, rewards are negative)
numVar = numEdges * 2 + numNodeConf + numRegions;
f = zeros(numVar,1); % init

edgeCostCoef = 100;
regionCostCoef = 500;            
%% Get initial training labels for the image

[labelImg_indexed,numLabels] = getLabelIndexImg(labelImage);
[c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,c_extNodeInds]...
            = getCells2WSregions(labelImg_indexed,ws,numLabels,setOfRegions,...
            edgeListInds,edges2nodes);

activeEdgeIDs = getElementsFromCell(c_extEdgeIDs);
[~,activeEdgeListInds] = intersect(edgeListInds,activeEdgeIDs);

activeWSregionListInds = getElementsFromCell(c_cells2WSregions);

%activeRegionListInds = activeWSregionListInds - 1;

activeNodeListInds = getElementsFromCell(c_extNodeInds);

% edge variables
edgeUnary = edgeUnary.*edgeCostCoef;
% reward both components of active (initial guess) edges
f(activeEdgeListInds) = -edgeUnary(activeEdgeListInds);
activeEdgeListInds_complementary = activeEdgeListInds + numEdges;
f(activeEdgeListInds_complementary) = -edgeUnary(activeEdgeListInds);
% penalize both components of inactive (initial guess) edges
edgeSeq = 1:numEdges;
inactiveEdgeListInds = setdiff(edgeSeq,activeEdgeListInds);
f(inactiveEdgeListInds) = -edgeUnary(inactiveEdgeListInds);
inactiveEdgeListInds_complementary = inactiveEdgeListInds + numEdges;
f(inactiveEdgeListInds_complementary) = -edgeUnary(inactiveEdgeListInds);

% region variables
offSet_regionVar = numVar - numRegions;

offSet_ActiveRegionLIDs = offSet_regionVar + activeWSregionListInds;
f(offSet_ActiveRegionLIDs) = -regionCostCoef; % reward

regionSeq = 1:numRegions;
inactiveWSregionListInds = setdiff(regionSeq,activeWSregionListInds);
offSet_InactiveWSregionLIDs = inactiveWSregionListInds + offSet_regionVar;
f(offSet_InactiveWSregionLIDs) = regionCostCoef; % penalty







