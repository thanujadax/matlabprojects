function neuronIDs2D = get2DComponentStructures(neuronIDsForGridCells,...
                        numR,numC,numZ)

% assign a unique ID for each 2D neuron section

% Outputs:
%   c_p2c: parentID -> childrenIDs
%   adjGraph: adjacency graph for all nodes(=2Dneuron section)
%           TODO: make this a sparse matrix
%   edgePropertices: edgeID -> properties

% Input:
%   neuronIDsForGridCells: N-by-1 matrix where each element says the
%   neuronID the gridCell belongs to.

% Init
[numR,numC,numZ] = size(neuronIDsForGridCells);
numGridCellsPerSlice = numR * numC;
numGridCells = numel(neuronIDsForGridCells);

maxNeuronID = max(max(max(neuronIDsForGridCells)));
NID2D = 0;

neuronIDs2D = zeros(numGridCells,1);
c_p2c = cell(maxNeuronID,1);
adjGraph = zeros(maxNeuronID,maxNeuronID);

gridStopInd = 0;

for k=2:(numZ-1)
    % get each 2D connected comp (by neuronID)
    gridStartInd = gridStopInd + 1;
    gridStopInd = gridStopInd + numGridCellsPerSlice;
    gridForThisSlice = neuronIDsForGridCells(gridStartInd:gridStopInd);
    slice_i = reshape(gridForThisSlice,numR,numC);
           
    % get what's connected in layer k+1 (if any)
    
    
    
end