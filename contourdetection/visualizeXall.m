function ilpSegmentation = visualizeXall(x,sizeR,sizeC,numEdges,numRegions,edgepixels,...
            junctionTypeListInds,nodeInds,connectedJunctionIDs,...
            nodeEdges,edgeListInds,wsIDsForRegions,ws,twoRegionEdges,edges2regions,...
            output,showIntermediate)

% visualize active labels of x without considering any structure.

% get active edges and active nodes from x
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

activeContourPixels = find(ilpSegmentation);

% figure;imagesc(ilpSegmentation);title('ILP contours');
% reconstruct the edges with the values from the orientation filters (HSV)
% [output rgbimg] = reconstructHSVgauss_mv(orientedScoreSpace3D,orientations,...
%             barLength,barWidth,threshFrac,medianFilterH);
% get the active pixels
%output(:,:,3) = ilpSegmentation;
totX = numel(x);
numRegionVars = numRegions * 2;
% get active foreground cells
regionStartPos = totX - numRegionVars + 2;
regionActivationVector = x(regionStartPos:2:totX);
activeRegionInd = find(regionActivationVector>0);
% get internal pixels for foreground cells
foregroundPixels = [];

% wsIDs = getWsIDsForCellIDs(ws,setOfCells,edges2pixels,nodeInds,...
%             edges2nodes,edgeListInds);

for i=1:numel(activeRegionInd)
    wsID_i = wsIDsForRegions(activeRegionInd(i));
    if(wsID_i==0)
        disp('problem with wsid check')
    else
%         cellPixInds_i = find(ws==wsID_i);
        cellPixInds_i = getInternalPixForCell(ws,wsID_i);
        foregroundPixels = [foregroundPixels; cellPixInds_i];
    end
end

% foregroundPixels = setdiff(foregroundPixels,offEdgePixelInds);

% make inactive edges inside foreground regions show as foreground
[inEdgeListInds,inEdgeIDs] = getInEdges(twoRegionEdges,regionActivationVector,...
                onEdgeStates,edges2regions,edgeListInds);
            
if(~isempty(inEdgeListInds))
    inEdgePixels = [];
    for i=1:numel(inEdgeListInds)
        inEdgePixels_i = edgepixels(inEdgeListInds(i),:);
        inEdgePixels_i = inEdgePixels_i(inEdgePixels_i>0);
        inEdgePixels = [inEdgePixels; inEdgePixels_i'];
    end
    foregroundPixels = [foregroundPixels; inEdgePixels];
end

% make inactive nodes inside foreground regions show as foreground
% nodeActivationVector = ~nodeInactiveStates_x;

inNodePixels = getInNodePixels(inEdgeIDs,nodeEdges,...
        nodeActivationVector,connectedJunctionIDs);
foregroundPixels = [foregroundPixels; inNodePixels];


foregroundPixels = setdiff(foregroundPixels,activeContourPixels);


% first visualize the contours
output_h = output(:,:,1);
output_s = output(:,:,2);
output_v = ilpSegmentation;

% create HSV image
% hsvImage = cat(3,output(:,:,1),output(:,:,2),ilpSegmentation);
hsvImage = cat(3,output_h,output_s,output_v);
% convert it to an RGB image
RGBimg = hsv2rgb(hsvImage);
% titleStr = sprintf('C = %d : lambda = %d',cNode,decayRate);
% titleStr = sprintf('C = %d',maxCost_direction);
if(showIntermediate)
    figure;imshow(RGBimg)
end

% assign white to active (foreground) cells

output_h(foregroundPixels) = 1; 
output_s(foregroundPixels) = 0;
output_v(foregroundPixels) = 1;

% create HSV image
% hsvImage = cat(3,output(:,:,1),output(:,:,2),ilpSegmentation);
hsvImage_foreground = cat(3,output_h,output_s,output_v);
% convert it to an RGB image
RGBimg_foreground = hsv2rgb(hsvImage_foreground);
% titleStr = sprintf('C = %d : lambda = %d',cNode,decayRate);
% titleStr = sprintf('C = %d',maxCost_direction);
if(showIntermediate)
    figure;imshow(RGBimg_foreground)
end
