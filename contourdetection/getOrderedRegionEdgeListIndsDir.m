function [c_edgeLIDsForRegions_dir_cw,setOfRegions_edgeLIDs,edgeLIDs2nodes_directional] ...
        = getOrderedRegionEdgeListIndsDir...
        (setOfRegions,edges2nodes,jAnglesAll_alpha,...
        junctionTypeListInds,nodeEdgeIDs,edgeListIndsAll,...
        edges2pixels,sizeR,sizeC)

% Inputs:
%   setOfRegions: edgeIDs for each region as a row vector

% Outputs:
%   c_edgeLIDsForRegions_cw: set of cells each containing directional
%   edgeLIDs for the corresponding region
%   setOfRegions_edgeLIDs: edgeLIDs for each region (undirected)

% at each node take a right turn. get the edge. check if that edge is in
% the set of edges

numRegions = size(setOfRegions,1);
c_edgeLIDsForRegions_dir_cw = cell(numRegions,1);
setOfRegions_edgeLIDs = setOfRegions;

edges2nodes_complements = edges2nodes;
edges2nodes_complements(:,1) = edges2nodes(:,2);
edges2nodes_complements(:,2) = edges2nodes(:,1);
edgeLIDs2nodes_directional = [edges2nodes; edges2nodes_complements];

for i=1:numRegions
    edgeIDsOfRegion_i = setOfRegions(i,:);
    edgeIDsOfRegion_i = edgeIDsOfRegion_i(edgeIDsOfRegion_i>0);
    % arrange the edges in clockwise order of node traversal
    [~,edgeLIDsForRegion] = intersect(edgeListIndsAll,edgeIDsOfRegion_i);
    if(numel(edgeLIDsForRegion)==0)
        c_edgeLIDsForRegions_dir_cw{i} = 0;
    else
        numE_region = numel(edgeLIDsForRegion);
        setOfRegions_edgeLIDs(i,1:numE_region) = edgeLIDsForRegion;

        nodeLIdsForRegion = edgeLIDs2nodes_directional(edgeLIDsForRegion,:);
        nodeLIdsForRegion = unique(nodeLIdsForRegion);
        c_edgeLIDsForRegions_dir_cw{i} = getCwOrderedEdgesForRegion...
                (edgeLIDsForRegion,edgeLIDs2nodes_directional,junctionTypeListInds,...
                nodeEdgeIDs,jAnglesAll_alpha,edgeListIndsAll,nodeLIdsForRegion,...
                edges2nodes,edges2pixels,sizeR,sizeC);
    end
    
end


function cwOrderedDirEdgeListInds = getCwOrderedEdgesForRegion...
            (edgeListInds_region,edges2nodes_directional,junctionTypeListInds,...
            nodeEdgeIDs,jAnglesAll_alpha,edgeListIndsAll,nodeLIdsForRegion,...
            edges2nodes,edges2pixels,sizeR,sizeC)
% Pick one edge (1st in the list)
edgeLId_1 = edgeListInds_region(1);
edgeID_1 = edgeListIndsAll(edgeLId_1);

% Get the nodes at each end of the edge. At each node get the next edge as
% if to complete a clockwise cycle.
nodeListInds = edges2nodes_directional(edgeLId_1,:);
nextCwEdgeLId_1 = getNextClockwiseEdge(nodeListInds(1),edgeLId_1,edgeID_1,...
            nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListIndsAll,...
            edges2pixels,sizeR,sizeC);

nextCwEdgeLId_2 = getNextClockwiseEdge(nodeListInds(2),edgeLId_1,edgeID_1,...
            nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListIndsAll,...
            edges2pixels,sizeR,sizeC);

% One of the two edges belongs to the current region. Keep this edge as the
% next edge
nextCwEdgeLId_inRegion = intersect...
                (edgeListInds_region,[nextCwEdgeLId_1,nextCwEdgeLId_2]);
            
            
% debug code start 20141224
if(numel(nextCwEdgeLId_inRegion)>1)
    disp('getOrderedRegionEdgeListIndsDir. numEdges >1.. running tie breaker...');
    nextCwEdgeLId_inRegion = nextCwEdgeTieBreaker(...
    nodeInd,prevEdgeLID,nextEdgeLIDsList,edgepixels,...
    connectedJunctionIDs,sizeR,sizeC);
    disp('tie broken');
    % error('ERROR: getOrderedRegionEdgeListIndsDir. numEdges >1')
elseif(numel(nextCwEdgeLId_inRegion)<1)
    error('ERROR: getOrderedRegionEdgeListIndsDir. numEdges <1')
    
else

    % debug code stop

    % % debug code start (???): important to make sure the outputs are assigned
    % i = 1;
    % % edgeListInds_region = edgeListInds_region(edgeListInds_region>0);
    % while(numel(nextCwEdgeLId_inRegion)~=1)
    %     % Pick one edge (1st in the list)
    %     
    %     if(numel(edgeListInds_region)>=i);
    %         edgeLId_1 = edgeListInds_region(i);
    %         edgeID_1 = edgeListIndsAll(edgeLId_1);
    %         % Get the nodes at each end of the edge. At each node get the next edge as
    %         % if to complete a clockwise cycle.
    %         nodeListInds = edges2nodes_directional(edgeLId_1,:);
    %         nextCwEdgeLId_1 = getNextClockwiseEdge(nodeListInds(1),edgeLId_1,edgeID_1,...
    %                     nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListIndsAll,...
    %                 edges2pixels,sizeR,sizeC);
    % 
    %         nextCwEdgeLId_2 = getNextClockwiseEdge(nodeListInds(2),edgeLId_1,edgeID_1,...
    %                     nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListIndsAll,...
    %                     edges2pixels,sizeR,sizeC);
    % 
    %         % One of the two edges belong to the current region. Keep this edge as the
    %         % next edge
    %         nextCwEdgeLId_inRegion = intersect...
    %                         (edgeListInds_region,[nextCwEdgeLId_1,nextCwEdgeLId_2]);   
    %     end
    %     i = i + 1;
    % end
    % % debug code end

    if(isempty(nextCwEdgeLId_inRegion))
        error('ERROR1: getOrderedRegionEdgeListIndsDir.m nextCwEdge not found!')
    else
        % From this edge, pick the node at the other end find the edge attached to
        % it in the same region. This is the next edge. Continue finding the next
        % edges in the clockwise order until all the edges are collected.
        if(nextCwEdgeLId_inRegion==nextCwEdgeLId_1)
            nodeLId_0 = nodeListInds(1);
        else
            nodeLId_0 = nodeListInds(2);
        end

        cwOrderedDirEdgeListInds = getCwDirSetOfEdges(edgeLId_1,nextCwEdgeLId_inRegion,...
                        nodeLId_0,edges2nodes_directional,edgeListInds_region,...
                        nodeLIdsForRegion,nodeEdgeIDs,edgeListIndsAll,edges2nodes,...
                        junctionTypeListInds,jAnglesAll_alpha,...
                        edges2pixels,sizeR,sizeC);

    end

end % new end added to match the new debug code 20141224

function cwOrderedDirEdgeListInds = getCwDirSetOfEdges(edgeLId_1,nextEdgeLID,...
                        nextNodeLID,edges2nodes_directional,edgeListInds_region,...
                        nodeLIdsForRegion,nodeEdgeIDs,edgeListIndsAll,edges2nodes,...
                        junctionTypeListInds,jAnglesAll_alpha,...
                        edges2pixels,sizeR,sizeC)       
                    
numEdges_region = numel(edgeListInds_region);
cwOrderedDirEdgeListInds = zeros(numEdges_region,1);
cwOrderedNodeListInds = zeros(numEdges_region,1);

% cwOrderedEdgeListInds(1) = edgeLId_1; % TODO: use directed LID
numEdges = numel(edgeListIndsAll);
nextEdgeNodeLID_pair = edges2nodes(nextEdgeLID,:);
if(nextEdgeNodeLID_pair(1)==nextNodeLID)
    nextEdgeLID_dir = nextEdgeLID;
else
    nextEdgeLID_dir = nextEdgeLID + numEdges;
end

cwOrderedDirEdgeListInds(1) = nextEdgeLID_dir; % do just (this) one

cwOrderedNodeListInds(1) = nextNodeLID;
nextNodePair = edges2nodes(nextEdgeLID,:);
N1 = nextNodeLID;
% N2 is at the other end of the edge edgeLId_2;
nodeLIds_temp = edges2nodes_directional(nextEdgeLID,:);
N2 = setdiff(nodeLIds_temp,N1);

%cwOrderedEdgeListInds(1) = edgeLId_2; % is it in the right direction
if(numEdges_region>1)
    for i = 2:numEdges_region
        % nextNodePair = edges2nodes(nextEdgeLID,:);
        nextNodeLID = setdiff(nextNodePair,nextNodeLID);
        cwOrderedNodeListInds(i) = nextNodeLID;

        allEdgeIDsForNextNode = nodeEdgeIDs(nextNodeLID,:);
        allEdgeIDsForNextNode(1) = []; % 1st element is nodePixInd
        [~,allEdgeLIDsForNextNode] = intersect(edgeListIndsAll,allEdgeIDsForNextNode);
        nodeEdgeLIDpair = intersect(edgeListInds_region,allEdgeLIDsForNextNode);
        if(numel(nodeEdgeLIDpair)>2)
            % depending on the incoming edge, pick the next node pair
            % get the ccw angles wrt incoming edge
            
            % get the edge with the closest ccw angle as the next edge
            nextEdgeLID = getNextClockwiseEdge(nextNodeLID,0,edgeListIndsAll(nextEdgeLID),...
            nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListIndsAll,...
            edges2pixels,sizeR,sizeC);
        
        elseif(numel(nodeEdgeLIDpair)<2)
           disp('ERROR: getOrderedRegionEdgeListIndsDir.m: no next edge!')  
            
        else            
            nextEdgeLID = setdiff(nodeEdgeLIDpair,nextEdgeLID); 
        end
        nextNodePair = edges2nodes(nextEdgeLID,:);     
        if(numel(nextNodePair)==0)
           disp('ERROR: getOrderedRegionEdgeListIndsDir.m: no next node pair!!') 
        end
        
        if(nextNodePair(1)==nextNodeLID)
            nextEdgeLID_dir = nextEdgeLID;
        else
            nextEdgeLID_dir = nextEdgeLID + numEdges;
        end
        if(numel(nextEdgeLID_dir) ~= 1)
            disp('ERROR: getOrderedRegionEdgeListIndsDir.m: numel(nextEdgeLID_dir) ~= 1');
            numel(nextEdgeLID_dir)
        else
        cwOrderedDirEdgeListInds(i) = nextEdgeLID_dir;
        end
    end
end

% for i=1:numEdges_region
%     % get all edges connected to N2 as the starting node (N1)
%     N1 = N2;
%     edgesForN1_logical = (edges2nodes_directional(:,1)==N1);
%     allN2sForN1 = edges2nodes_directional(edgesForN1_logical,2);
%     % which one of those edges belongs to this region
%     N2 = intersect(allN2sForN1,nodeLIdsForRegion);
%     % which edge is N1->N2: which row has (N1,N2)?
%     nextEdgeLID_dir = strmatch([N1,N2],edges2nodes_directional);
%     cwOrderedEdgeListInds(i) = nextEdgeLID_dir;    
% end


function nextCwEdgeLInd = getNextClockwiseEdge(nodeLId,edgeLId,edgeID,...
            nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListInds,...
            edges2pixels,sizeR,sizeC)
% at the node nodeLId, wrt the edge edgeLId, what is the next edge in cw direction
nodeEdgeIDsAll = nodeEdgeIDs(nodeLId,:);
nodeEdgeIDsAll(1) = []; % first element is the nodePixelInd
nodeEdgeIDsAll = nodeEdgeIDsAll(nodeEdgeIDsAll>0);

[junctionListInd,junctionType] = find(junctionTypeListInds==nodeLId);
alphas_junctionType = jAnglesAll_alpha{junctionType};

nodeAlphas_0 = alphas_junctionType(junctionListInd,:);
numOfUniqueAngles = unique(nodeAlphas_0);
if(numel(numOfUniqueAngles)<numel(nodeAlphas_0))
    % duplicate alphas detected. recalculate using only the immediate edge pixel
    % wrt to node pixel
    nodeAlphas_0 = recalculateAlphas(nodeLId,nodeEdgeIDs,...
        edges2pixels,sizeR,sizeC);
end
% inputEdgePos = find(nodeEdgesAll==edgeID);
alpha_inputEdge = nodeAlphas_0(nodeEdgeIDsAll==edgeID);
% clockwise search
% get the edge with next smaller alpha
% get all angles smaller than alpha_0    
smallerAlphas = nodeAlphas_0(nodeAlphas_0<alpha_inputEdge);
if(~isempty(smallerAlphas))
    nextAlpha = max(smallerAlphas);
else
    % get the largest alpha
    nextAlpha = max(nodeAlphas_0);
end
nextCwEdgeID = nodeEdgeIDsAll(nodeAlphas_0==nextAlpha);
[~,nextCwEdgeLInd] = intersect(edgeListInds,nextCwEdgeID);


function alphas = recalculateAlphas(nodeListInd,nodeEdges,edges2pixels,sizeR,sizeC)
% alphas wrt to the node pixel using only the immediate edge pixels
% inputs:
%   nodeListInd
%   nodeEdges
%   edges2pixels
MAX_NUM_PIXELS = 1;
nodePixInd = nodeEdges(nodeListInd,1);
nodeEdgeIDs = nodeEdges(nodeListInd,:);
nodeEdgeIDs(1) = []; % the first element is node
nodeEdgeIDs = nodeEdgeIDs(nodeEdgeIDs>0);
numEdges = numel(nodeEdgeIDs);
alphas = zeros(1,numEdges);
[y0,x0] = ind2sub([sizeR,sizeC],nodePixInd);

for i=1:numEdges
    edgePixelInds=edges2pixels((edges2pixels(:,1)==nodeEdgeIDs(i)),:);
    edgePixelInds(1) = []; % first element is the edgeID
    edgePixelInds = edgePixelInds(edgePixelInds>0);
    nodeEdgePixel = getNodeEdgePixel(nodePixInd,edgePixelInds,sizeR,sizeC,...
                                    MAX_NUM_PIXELS);
    [y1,x1] = ind2sub([sizeR,sizeC],nodeEdgePixel);
    y = y1 - y0;
    x = x1 - x0;
    
    alpha_i = atan2d(y,x);
    if(alpha_i<0)
        alpha_i = alpha_i + 360;
    end
    
    alphas(i) = alpha_i;

end % for loop
% end of function recalculateAlphas