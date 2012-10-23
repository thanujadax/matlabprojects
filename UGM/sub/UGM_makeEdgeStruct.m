function [edgeStruct] = UGM_makeEdgeStruct(adj,nStates,useMex,maxIter)
% [edgeStruct] = UGM_makeEdgeStruct(adj,nStates,useMex,maxIter)
%
% adj - nNodes by nNodes adjacency matrix (0 along diagonal)
%

if nargin < 3
    useMex = 1;
end
if nargin < 4
    maxIter = 100;
end

nNodes = int32(length(adj));
[i j] = ind2sub([nNodes nNodes],find(adj));
nEdges = length(i)/2;
edgeEnds = zeros(nEdges,2,'int32');
eNum = 0;
for e = 1:nEdges*2
   if j(e) < i(e)
       edgeEnds(eNum+1,:) = [j(e) i(e)];
       eNum = eNum+1;
   end
end

[V,E] = UGM_makeEdgeVE(edgeEnds,nNodes,useMex);
% E - [ [indexes of edges connected to node 1][indexes of edges connected to node 2]
% .... [indexes of edges connected to node nNodes] ]
% V -- V(i) = the sum of the number of edges connected to nodes (1,2,...i-1)
% plus 1  (the sums of the lengths of the above blocks 1,2..i-1 plus 1)

edgeStruct.edgeEnds = edgeEnds;
edgeStruct.V = V;
edgeStruct.E = E;
edgeStruct.nNodes = nNodes;
edgeStruct.nEdges = size(edgeEnds,1);

% Handle other arguments
if isscalar(nStates)
   nStates = repmat(nStates,[double(nNodes) 1]);
end
edgeStruct.nStates = int32(nStates(:));
edgeStruct.useMex = useMex;
edgeStruct.maxIter = int32(maxIter);


