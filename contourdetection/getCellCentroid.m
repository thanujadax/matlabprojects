function cellCentroid = getCellCentroid(edgeIDset_cell,edges2pixels,edgeIDs_all,...
    sizeR,sizeC,edges2nodes,nodeInds)

% Inputs:
%   wsgraph - matrix containing non-zero valued pixels for edges and
%   junctions. zeros otherwise
%   edgeIDset - set of edges defining the current cell
%   edges2pixels - first column shows the actual edgeID. each row has the
%   list of pixels for the corresponding edge

% Outputs:
%   cellCentroid = [x; y] coordinates of the approximate centroid of the
%   cell

%% method 1
% % get all the edges bounding the cell
% numEdgesForCell = numel(edgeIDset_cell);
% edgeListInds = zeros(numEdgesForCell,1);
% X = [];
% Y = [];
% for i=1:numEdgesForCell
%     edgeListInds(i) = find(edgeIDs_all==edgeIDset_cell(i));
%     clear edge_i_pixels
%     edge_i_pixels = edges2pixels(edgeListInds(i),:);
%     edge_i_pixels(1) = []; % first element is the edgeID
%     edge_i_pixels = edge_i_pixels(edge_i_pixels>0);
%     clear x y
%     [y x] = ind2sub([sizeR sizeC],edge_i_pixels');
%     X = [X; x];
%     Y = [Y; y];
% end
% 
% meanx = floor(mean(X));
% meany = floor(mean(Y));
% 
% cellCentroid = [meanx meany];

%% method 2
% get the vertices (x,y)
numEdgesForCell = numel(edgeIDset_cell);
edgeListInds = zeros(numEdgesForCell,1);
nodeset = [];
for i=1:numEdgesForCell
    edgeListInds(i) = find(edgeIDs_all==edgeIDset_cell(i));
    nodes = edges2nodes(edgeListInds(i),:);
    nodeset = [nodeset nodes];
end
nodeset = unique(nodeset);
nodeset_pixinds = nodeInds(nodeset);
[y,x] = ind2sub([sizeR sizeC],nodeset_pixinds');
N=numel(x);
y(end+1) = y(1);
x(end+1) = x(1);
% calculate area A
sum = 0;
for i=1:N
    sum = sum + x(i)*y(i+1) - x(i+1)*y(i);
end
A = sum/2;
% cog
sum = 0;
for i=1:N
    sum = sum + (x(i) + x(i+1)) * (x(i)*y(i+1) - x(i+1)*y(i));
end
X = sum/(6*A);
sum = 0;
for i=1:N
    sum = sum + (y(i) + y(i+1)) * (x(i)*y(i+1) - x(i+1)*y(i));
end
Y = sum/(6*A);

cellCentroid = [floor(X) floor(Y)];