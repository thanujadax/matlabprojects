function adjacencyMat = getNodeAdjacency_sq(nodePix)

% returns adjacency matrix for the given nodePix laid out on a triagular
% grid

% input:
%    nodePix - sizeR x sizeC matrix. Conserves the order of the square grid
%    layout (adjacencies are conserved).

% get boundary nodes

% for non-boundary nodes
% get neighbors at (x,y+d), (x,y-d), (x+d,y), (x-d,y)


numNodes = numel(nodePix);
nodeIndsVect = reshape(nodePix,numNodes,1);
% TODO: make adjacencyMat a sparse matrix
adjacencyMat = zeros(numNodes);

% for all non boundary grid nodes
% row updates of adjacency matrix
for i=2:numNodes
    for j=2:numNodes
        nodePix_i = nodePix(i,j);
        % current pix
        nodeListInd_i = find(nodeIndsVect==nodePix_i);
        % get 4-neighbors
        fourN = get4N(i,j,nodePix);
        % update row for nodeListInd_i
        adjacencyMat(nodeListInd_i,fourN) = nodeListInd_i;
    end
end

%% for boundary grid nodes (except the corners)
% 1 - top row
i=1;
for j=2:(sizeC-1)
    isTopRow = 1;
    threeN = get3N_LRD(i,j,nodePix,isTopRow);
    % get NodeListInd_i (row ID to be updated)
    nodePix_i = nodePix(i,j);
    % current pix
    nodeListInd_i = find(nodeIndsVect==nodePix_i);    
    % update adjacencyMat row
    adjacencyMat(nodeListInd_i,threeN) = nodeListInd_i;
end
% 2 - bottom row
i=sizeR;
for j=2:(sizeC-1)
    isTopRow = 0;
    threeN = get3N_LRD(i,j,nodePix,isTopRow);
    % get NodeListInd_i (row ID to be updated)
    nodePix_i = nodePix(i,j);
    % current pix
    nodeListInd_i = find(nodeIndsVect==nodePix_i);    
    % update adjacencyMat row
    adjacencyMat(nodeListInd_i,threeN) = nodeListInd_i;
end
% 3 - left column
j=1;
for i=2:(sizeR-1)
    isRight = 0;
    threeN = get3N_UDR(i,j,nodePix,isRight);    
    % get NodeListInd_i (row ID to be updated)
    nodePix_i = nodePix(i,j);
    % current pix
    nodeListInd_i = find(nodeIndsVect==nodePix_i);    
    % update adjacencyMat row
    adjacencyMat(nodeListInd_i,threeN) = nodeListInd_i;
end
% 4 - right column
j=sizeC;
for i=2:(sizeR-1)
    isRight = 1;
    threeN = get3N_UDR(i,j,nodePix,isRight);    
    % get NodeListInd_i (row ID to be updated)
    nodePix_i = nodePix(i,j);
    % current pix
    nodeListInd_i = find(nodeIndsVect==nodePix_i);    
    % update adjacencyMat row
    adjacencyMat(nodeListInd_i,threeN) = nodeListInd_i;
end

%%  for the four corners
% 1 - top left corner
twoN = zeros(2,1);
i = 1;
j=1;
x = i;
y = j+1;
twoN(1) = nodePix(x,y);
x = i+1;
y = j;
twoN(2) = nodePix(x,y);
nodePix_i = nodePix(i,j);
% current pix
nodeListInd_i = find(nodeIndsVect==nodePix_i);    
% update adjacencyMat row
adjacencyMat(nodeListInd_i,twoN) = nodeListInd_i;

% 2 - bottom left corner
twoN = zeros(2,1);
i = sizeR;
j=1;
x = i;
y = j+1;
twoN(1) = nodePix(x,y);
x = i-1;
y = j;
twoN(2) = nodePix(x,y);
nodePix_i = nodePix(i,j);
% current pix
nodeListInd_i = find(nodeIndsVect==nodePix_i);    
% update adjacencyMat row
adjacencyMat(nodeListInd_i,twoN) = nodeListInd_i;

% 3 - bottom right corner
twoN = zeros(2,1);
i = sizeR;
j=1;
x = i;
y = j-1;
twoN(1) = nodePix(x,y);
x = i-1;
y = j;
twoN(2) = nodePix(x,y);
nodePix_i = nodePix(i,j);
% current pix
nodeListInd_i = find(nodeIndsVect==nodePix_i);    
% update adjacencyMat row
adjacencyMat(nodeListInd_i,twoN) = nodeListInd_i;

% 4 - top right corner
twoN = zeros(2,1);
i = sizeR;
j=1;
x = i;
y = j-1;
twoN(1) = nodePix(x,y);
x = i+1;
y = j;
twoN(2) = nodePix(x,y);
nodePix_i = nodePix(i,j);
% current pix
nodeListInd_i = find(nodeIndsVect==nodePix_i);    
% update adjacencyMat row
adjacencyMat(nodeListInd_i,twoN) = nodeListInd_i;

%% supplementary functions
function threeN = get3N_LRD(i,j,nodePix,isTopRow)
% left right 2N and below/above
threeN = zeros(2,1);
x = i;
y = j-1;
threeN(1) = nodePix(x,y);

x = i;
y = j+1;
threeN(2) = nodePix(x,y);

if(isTopRow)
x = i+1;
else
x = i-1;    
end
y=j;
threeN(3) = nodePix(x,y);


function threeN = get3N_UDR(i,j,nodePix,isRight)
% up down 2N and left/right
threeN = zeros(2,1);
x = i-1;
y = j;
threeN(1) = nodePix(x,y);

x = i+1;
y = j;
threeN(2) = nodePix(x,y);

if(isRight)
    y = j-1;
else
    y = j + 1;
end
x = i;
threeN(3) = nodePix(x,y);

function fourN = get4N(i,j,nodePix)    
    fourN = zeros(4,1);
    x = i;
    y = j-1;
    fourN(1) = nodePix(x,y);
    
    x = i;
    y = j+1;
    fourN(2) = nodePix(x,y);
    
    x = i-1;
    y = j;
    fourN(3) = nodePix(x,y);
    
    x = i+1;
    y = j;
    fourN(4) = nodePix(x,y);
    