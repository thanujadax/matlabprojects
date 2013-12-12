function [c_wsIDsInCell,c_internalEdgeLIDsInCell,c_extDirectedEdgeLIDsInCell,...
          c_internalNodeListInds,c_extNodeListInds]...
            = getCells2WSregions2(labelImg_indexed,ws,numLabels,setOfRegions,...
            edgeListInds,edges2nodes,junctionTypeListInds,nodeInds,edgepixels)

% version 2. with directed edges. outputs edgeLIDs instead of eIDs.
        
% Inputs:
%   setOfRegions: contains edgeIDs for each wsRegion
%   nodeEdges: edgeIDs for each nodePixelInd
%   connectedJunctionIDs: 
%   edges2nodes: gives 2 nodeLIDs for each edgeID

% Outputs:
%     c_wsIDsInCell
%     c_internalEdgeIDsInCell
%     c_extEdgeIDsInCell
%     c_internalNodeListInds
%     c_extNodeListInds


% Params:
threshFrac = 0.33; % fraction of WSregion pixels allowed to be OUTSIDE of the cell

c_wsIDsInCell = cell(numLabels,1);
c_internalEdgeLIDsInCell = cell(numLabels,1); % off edges - both directions
c_extDirectedEdgeLIDsInCell = cell(numLabels,1); % boundary edges LIDs
c_internalNodeListInds = cell(numLabels,1); % off nodes
c_extNodeListInds = cell(numLabels,1);

% regionID = wsID - 1;


for i=1:numLabels
    % label_i corresponds to cell_i
    ws_i = ws;
    setOfRegions_i = setOfRegions;
    pixelsForLabel_i_logical = (labelImg_indexed==i);
    wsIDsForLabel_i = ws_i(pixelsForLabel_i_logical);
    wsIDsForLabel_i = unique(wsIDsForLabel_i);
    wsIDsForLabel_i = wsIDsForLabel_i(wsIDsForLabel_i>0);
    
    % pick the wsRegions having 4/5 of their pixels in this cell
    wsIDsForLabel_i_filtered = getGoodWSregions...
        (wsIDsForLabel_i,ws,pixelsForLabel_i_logical,threshFrac);
    
    c_wsIDsInCell{i} = wsIDsForLabel_i_filtered; 
    regionIDs = wsIDsForLabel_i_filtered - 1;
    regionIDs = regionIDs(regionIDs>0);
    % get all edges for the cell
    edgesForCell_i = setOfRegions_i(regionIDs,:);
    edgesForCell_i = edgesForCell_i(edgesForCell_i>0);
    
    
    % pick the internal ones, having 2 regions bounded by it, which are
    % part of this cell
    edgeIDs_unique_i = unique(edgesForCell_i);
    edgeCounts_i = histc(edgesForCell_i,edgeIDs_unique_i);
    internalEdgeIDs_i = edgeIDs_unique_i(edgeCounts_i>1);
    [~,internalEdgeLIDs_i] = intersect(edgeListInds,internalEdgeIDs_i);
    c_internalEdgeLIDsInCell{i} = internalEdgeLIDs_i;
    
    edges2nodes_complements = edges2nodes;
    edges2nodes_complements(:,1) = edges2nodes(:,2);
    edges2nodes_complements(:,2) = edges2nodes(:,1);
    edges2nodes_directed = [edges2nodes; edges2nodes_complements];
    
    % the other edges are external
    extEdgeIDs_i = setdiff(edgeIDs_unique_i,internalEdgeIDs_i);
    
    % all nodes where 2 internal edges meet are internal nodes    
    % nodes where at least 1 external edge meets with other nodes, are
    % external nodes
    [internalNodeListInds,extNodeListInds] = getNodes...
                (internalEdgeIDs_i,extEdgeIDs_i,edges2nodes,edgeListInds);
    
    [~,extEdgeLIDs_i] = intersect(edgeListInds,extEdgeIDs_i);
    [cwOrderedEdgeLIDs,cwOrderedNodeLIDs] = getOrderedEdgeLIDsCw(extNodeListInds,...
            edges2nodes_directed,edgeListInds,extEdgeLIDs_i,...
            edges2nodes,junctionTypeListInds);
    
    c_extDirectedEdgeLIDsInCell{i} = cwOrderedEdgeLIDs; 
    
    
    c_internalNodeListInds{i} = internalNodeListInds;
    
    c_extNodeListInds{i} = cwOrderedNodeLIDs;
                
end



function [cwOrderedEdgeLIDs,cwOrderedNodeLIDs] = getOrderedEdgeLIDsCw(nodeList_cell,...
            edges2nodes_directed,edgeListInds,edgeListInds_cell,...
            edges2nodes,junctionTypeListInds)
% order the nodes and edges in cw order
% outputs cw_directed_edge_LIDs

% pick up first edgeLID_dir_cell
nextEdgeLID = edgeListInds_cell(1); 
nodesOfNextEdge = edges2nodes_directed(nextEdgeLID,:);
% pick up one node of this edge
numNodes_cell = numel(nodeList_cell);
numEdges_cell = numNodes_cell;

% to store ordered edgeLIDs and nodeLIDs
orderedEdgeLIDs_cell = zeros(numEdges_cell,1); % initialize
orderedNodeLIDs_cell = zeros(numNodes_cell,1);

orderedNodeLIDs_cell(1) = nodesOfNextEdge(1);
nextNodeLID = nodesOfNextEdge(2);
orderedNodeLIDs_cell(2) = nextNodeLID;

orderedEdgeLIDs_cell(1) = nextEdgeLID;

for i=2:numEdges_cell    
    % get the two edges connected to the other node
    [connectedEdgeIDs_toCurrentNode,~] = find(edges2nodes==nextNodeLID);
    % nextEdge is the other connected to the current node    
    [~,connectedEdgeLIDs_toCurrentNode] = intersect...
                        (edgeListInds,connectedEdgeIDs_toCurrentNode);
    % keep the one that is part of this cell
    connectedEdgeLIDs_toCurrentNode = intersect...
                    (edgeListInds_cell,connectedEdgeLIDs_toCurrentNode);
                
    nextEdgeLID = setdiff(connectedEdgeLIDs_toCurrentNode,nextEdgeLID);
    nextEdgeID = edgeListInds(nextEdgeLID);
    % nextNode is the other node of the next edge
    nextNodes_tmp = edges2nodes(nextEdgeID,:);
    nextNodeLID = setdiff(nextNodes_tmp,nextNodeLID);
    orderedNodeLIDs_cell(i+1) = nextNodeLID;
    
    orderedEdgeLIDs_cell(i) = nextEdgeLID;
end

% determine if the order is cw (or ccw)
isCw = checkEdgeLIDsCw(orderedNodeLIDs_cell,orderedEdgeLIDs_cell,...
            junctionTypeListInds,jEdges,jAnglesAll_alpha,edgeListInds);

% reverse the order if ccw
[cwOrderedEdgeLIDs,cwOrderedNodeLIDs,directedNodeListN1N2] = orderEdgesNodes...
            (isCw,orderedNodeLIDs_cell,edges2nodes_directed);


function [cwOrderedEdgeLIDs,cwOrderedNodeLIDs,directedNodeListN1N2]...
            = orderEdgesNodes...
            (isCw,orderedNodeLIDs_cell,edges2nodes_directed)

if(~isCw)
    % is ccw. reverse order
     % orderedEdgeLIDs_cell = flipud(orderedEdgeLIDs_cell);
     orderedNodeLIDs_cell = flipud(orderedNodeLIDs_cell);
end
[cwOrderedEdgeLIDs,cwOrderedNodeLIDs,directedNodeListN1N2] = getDirectionalEdgeLIDs...
        (orderedNodeLIDs_cell,edges2nodes_directed);


function [directedEdgeLIDs, directedNodeListN1N2] = getDirectionalEdgeLIDs...
        (orderedNodeLIDs_cell,edges2nodes_directional)
% append a copy of the first node at the end of the list
orderedNodeLIDs_cell(end+1) = orderedNodeLIDs_cell(1); 
numNodes_cell = numel(orderedNodeLIDs_cell);
numEdges_cell = numNodes_cell - 1;

directedEdgeLIDs = zeros(numEdges_cell,1);
directedNodeListN1N2 = zeros(numEdges_cell,2);

for i=1:numEdges_cell
    directedNodeListN1N2(i,1) = orderedNodeLIDs_cell(i);
    directedNodeListN1N2(i,2) = orderedNodeLIDs_cell(i+1);
    
    N1_logical = (edges2nodes_directional==directedNodeListN1N2(i,1));
    N2_logical = (edges2nodes_directional==directedNodeListN1N2(i,2));
    
    sum_logical_ind = N1_logical + N2_logical;
    
    directedEdgeLIDs(i) = find(sum_logical_ind==2);    
end

      
        
function isCw = checkEdgeLIDsCw(orderedNodeLIDs_cell,orderedEdgeLIDs_cell,...
            junctionTypeListInds,jEdges,jAnglesAll_alpha,edgeListInds)

numNodes_cell = numel(orderedNodeLIDs_cell);
increaseOfDirection = zeros(numNodes_cell);
orderedEdgeIDs_cell = edgeListInds(orderedEdgeLIDs_cell);

% rearrange the set of nodes
% N1 should go to N_end
nodeList_tmp = orderedNodeLIDs_cell;
nodeList_tmp(1:(end-1)) = orderedNodeLIDs_cell(2:end);
nodeList_tmp(end) = orderedNodeLIDs_cell(1);
orderedNodeLIDs_cell = nodeList_tmp;

for i=1:numNodesCell
    [nodeTypeLID,jType] = find(junctionTypeListInds==orderedNodeLIDs_cell(i));
    
    % get edgeIn and edgeOut
    edgeID_in = orderedEdgeIDs_cell(i);
    if(i<numNodesCell)
        edgeID_out = orderedEdgeIDs_cell(i+1);
    else
        edgeID_out = orderedEdgeIDs_cell(1);
    end
    
    nodeEdgeIDsAll_jType = jEdges{jType};
    nodeEdgeIDsAll_i = nodeEdgeIDs(nodeTypeLID,:);  
    
    alpha_jType = jAnglesAll_alpha{jType};
    alphasAll_node = alpha_jType(nodeTypeLID);
    
    % get the position of edgeIn and edgeOut in nodeEdges vector
    pos_edgeIn_logical = (nodeEdgeIDsAll_i==edgeID_in);
    pos_edgeOut_logical = (nodeEdgeIDsAll_i==edgeID_out);
    
    % get alpha1(edgeIn) and alpha2(edgeOut)
    alphaIn = alphasAll_node(pos_edgeIn_logical);
    % gammaIn is the leaving angle of edgeIn at the current node
    gammaIn = getGamma(alphaIn);
    alphaOut = alphasAll_node(pos_edgeOut_logical);
    
    % get the increase of direction and store in array
    increaseOfDirection(i) = alphaOut - gammaIn;
end

numIncreases = sum(increaseOfDirection>0);

fractionIncreases = numIncreases/numNodes_cell;
if(fractionIncreases>=0)
    isCw = 1;
else
    isCw = 0;
end


function gammaAngle = getGamma(alphaAngle)
if(alphaAngle >= 180)
    gammaAngle = alphaAngle - 180;
else
    gammaAngle = alphaAngle + 180;
end

function wsIDs_filtered = getGoodWSregions(wsIDsForLabel_i,ws,...
                            pixForLabel_j_logical,threshFrac)
wsIDs_filtered = [];
numWsRegions = numel(wsIDsForLabel_i);
for i=1:numWsRegions
    wsID_i = wsIDsForLabel_i(i);
    wsPixels_i_logical = (ws==wsID_i);
    wsPixOutsideLabelj_logical = ((pixForLabel_j_logical - wsPixels_i_logical)==-1);
    numWsPixOutside = sum(sum(wsPixOutsideLabelj_logical));
    numWsPixTot = sum(sum(wsPixels_i_logical));
    fractionOutside = numWsPixOutside/numWsPixTot;
    if (fractionOutside<threshFrac)
        wsIDs_filtered = [wsIDs_filtered; wsID_i];
    end
end

function [internalNodeListInds,extNodeListInds] = getNodes...
                (internalEdgeIDs,extEdgeIDs,edges2nodes,edgeListInds)

intEdgeListInd_logical = ismember(edgeListInds,internalEdgeIDs);            
internalNodeListInds_all = edges2nodes(intEdgeListInd_logical,:);
internalNodeListInds_all = unique(internalNodeListInds_all);

extEdgeListInd_logical = ismember(edgeListInds,extEdgeIDs);            
extNodeListInds = edges2nodes(extEdgeListInd_logical,:);
extNodeListInds = unique(extNodeListInds);

internalNodeListInds = setdiff(internalNodeListInds_all,extNodeListInds);

