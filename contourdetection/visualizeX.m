function segmentationOut = visualizeX(x,sizeR,sizeC,numEdges,numRegions,edgepixels,...
            junctionTypeListInds,nodeInds,connectedJunctionIDs,edges2nodes,...
            nodeEdges,edgeListInds,faceAdj,setOfRegions,wsIDsForRegions,ws,...
            marginSize,showIntermediate)

ilpSegmentation = zeros(sizeR,sizeC);
% active edges
% consider the edgeID given in the first col of edges2pixels?? no need for
% this since we use edgepixels array which is already sans the skipped

% edges
onStateEdgeXind = 2:2:(numEdges*2);
onEdgeStates = x(onStateEdgeXind);
onEdgeInd = find(onEdgeStates>0.5);
offEdgeListInd = find(onEdgeStates<0.5);
onEdgePixelInds = getPixSetFromEdgeIDset(onEdgeInd,edgepixels);
offEdgePixelInds = getPixSetFromEdgeIDset(offEdgeListInd,edgepixels);
ilpSegmentation(onEdgePixelInds) = 1;

% active nodes 
fIndStop = 2*numEdges;
nodeInactiveStates_x = [];
nodeActivationVector = zeros(numel(nodeInds),1);    % stores 1 for active node list inds
nodeIndsActive = [];
numJtypes = size(junctionTypeListInds,2);
for i=1:numJtypes
    % for each junction type
    % get the list of junctions and check their states in vector 'x'
    junctionNodesListListInds_i = find(junctionTypeListInds(:,i));
    if(~isempty(junctionNodesListListInds_i))
        junctionNodesListInds_i = junctionTypeListInds(junctionNodesListListInds_i,i);
        numJnodes_i = numel(junctionNodesListInds_i);
        % get the indices (wrt x) for the inactivation of the junctions
        numEdgePJ_i = i+1;
        numStatePJ_i = nchoosek(numEdgePJ_i,2)+1; % 1 is for the inactive case
        fIndStart = fIndStop + 1;
        fIndStop = fIndStart -1 + numJnodes_i*numStatePJ_i;
        fIndsToLook = fIndStart:numStatePJ_i:fIndStop; % indices of inactive state
        inactiveness_nodes_i = x(fIndsToLook);
        nodeInactiveStates_x = [nodeInactiveStates_x; inactiveness_nodes_i];
        activeStateNodeListInd = find(inactiveness_nodes_i<0.5);
        if(~isempty(activeStateNodeListInd))
            nodeListInd_i = junctionNodesListInds_i(activeStateNodeListInd);
            nodeActivationVector(nodeListInd_i) = 1;
            nodeIndsActive_i = nodeInds(nodeListInd_i);
            % if any of the active nodes are in the connectionJunction set,
            % make the other nodes in the same set active as well.
            for j=1:numel(nodeIndsActive_i)
                if(~isempty(connectedJunctionIDs))
                    indx = find(connectedJunctionIDs(:,1)==nodeIndsActive_i(j));
                    if(~isempty(indx))
                        % this is one of the cluster pixels
                        clusLabel = connectedJunctionIDs(indx,2);
                        clusNodeListInds = find(connectedJunctionIDs(:,2)==clusLabel); 
                        clusNodes = connectedJunctionIDs(clusNodeListInds,1);
                        ilpSegmentation(clusNodes) = 1;
                    end
                end
                
            end
            ilpSegmentation(nodeIndsActive_i) = 1;
            nodeIndsActive = [nodeIndsActive; nodeIndsActive_i];
        end
    end
end

inactiveNodePixInds = setdiff(nodeInds,nodeIndsActive);
offNodeIndList = find(ismember(nodeInds,inactiveNodePixInds));

totX = numel(x);
numRegionVars = numRegions*2;
% get active foreground cells
regionStartPos = totX - numRegionVars + 2;
regionActivationVector = x(regionStartPos:2:totX);
activeRegionInd = find(regionActivationVector>0);

%% store extracted geometry in datastructures
[c_cellBorderEdgeIDs,c_cellBorderNodeIDs] = getCellBorderComponents(onEdgeInd,...
        edges2nodes,nodeEdges,edgeListInds);

cellBorderPixels = getCellBorderPixels(c_cellBorderEdgeIDs,...
            c_cellBorderNodeIDs,edgepixels,nodeInds,connectedJunctionIDs);
        
visualizeCellBorders = zeros(sizeR,sizeC);
visualizeCellBorders(cellBorderPixels) = 1;
if(showIntermediate)
    figure;imshow(visualizeCellBorders)
end
% regions aggregating to form cells
% offEdgeIDList = edgeListInds(ismember(edgeListInds,offEdgeListInd)); 
offEdgeIDList = edgeListInds(offEdgeListInd);

[c_cells2regions,c_cellInternalEdgeIDs,c_cellIntNodeListInds] = getRegionsForOnCells(...
            faceAdj,activeRegionInd,offEdgeIDList,setOfRegions,wsIDsForRegions,ws,...
            offNodeIndList,edges2nodes,edgeListInds);
                
% visualize each cell in different colors
visualizeCells = zeros(sizeR,sizeC,3);
numCs = numel(c_cells2regions);
rMat = zeros(sizeR,sizeC);
gMat = zeros(sizeR,sizeC);
bMat = zeros(sizeR,sizeC);
for i=1:numCs
    % pick random color (RGB vals)
    R = rand(1); G = rand(1); B = rand(1);
    
    % get regions for this cell and the internal pixels. set RGB
    cellRegionList_i = c_cells2regions{i};
    regionPixels = getRegionPixels(cellRegionList_i,wsIDsForRegions,ws);
    rMat(regionPixels) = R;
    gMat(regionPixels) = G;
    bMat(regionPixels) = B;
    
    
    % get internal edges for this cell. set RGB.
    regionIntEdgeIDs = c_cellInternalEdgeIDs{i};
    regionIntEdgeListInds_logicalInd = ismember(edgeListInds,regionIntEdgeIDs);
    regionIntEdgePixels = edgepixels(regionIntEdgeListInds_logicalInd,:);
    regionIntEdgePixels = regionIntEdgePixels(regionIntEdgePixels>0);
    
    rMat(regionIntEdgePixels) = R;
    gMat(regionIntEdgePixels) = G;
    bMat(regionIntEdgePixels) = B;
      
    % get internal nodes for this cell. set RGB
    regionIntNodeListInds = c_cellIntNodeListInds{i};
    regionIntNodePixInds = getNodePixelsFromNodeInd...
                (regionIntNodeListInds,nodeInds,connectedJunctionIDs);
    
    rMat(regionIntNodePixInds) = R;
    gMat(regionIntNodePixInds) = G;
    bMat(regionIntNodePixInds) = B;
    
    visualizeCells(:,:,1) = rMat;
    visualizeCells(:,:,2) = gMat;
    visualizeCells(:,:,3) = bMat;
            
    
end

% figure;imshow(visualizeCells);
segmentationOut = removeThickBorder(visualizeCells,marginSize);

figure;imshow(segmentationOut);