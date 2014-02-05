function edges2pixels = getEdges2pixels_grid(edges2nodes,nodeIndsVect,...
                sizeR,sizeC,maxEdgeLenPix)

% first column contains the edgeID (=edgeLID in the non-grid versions)
            
% for each edge, linearly interpolate between the two nodes to get the
% intermediate pixels (=edge pixels)

numEdges = size(edges2nodes,1);

edges2pixels = zeros(numEdges,(maxEdgeLenPix+1));

for i=1:numEdges
    % get the two end points (nodes) for this edge
    nodes_i = edges2nodes(i,:);
    nodeInds_i = nodeIndsVect(nodes_i);
    % interpolate linearly to get the edge pixels between the nodes
    % intermediate points to evaluate: wrt x or y
    [y,x] = ind2sub([sizeR sizeC],nodeInds_i);
    yd = abs(y(1)-y(2));
    xd = abs(x(1)-x(2));
    if(yd>xd)
        % define reference points in r (y)
        yq = (min(y)+1) : (max(y)-1);
        xq = interp1(y,x,yq); % linear interpolation
    else
        % define reference points in c (x)
        xq = (min(x)+1) : (max(x)-1);
        yq = interp1(x,y,xq); % linear interpolation
    end
    % update edges2pixels
    edgePix_i = sub2ind([sizeR sizeC],yq,xq);
    edgeLen_i = numel(edgePix_i);
    edges2pixels(i,1) = i; % edgeID (=edgeLID)
    edges2pixels(i,2:(edgeLen_i+1)) = edgePix_i; 
end