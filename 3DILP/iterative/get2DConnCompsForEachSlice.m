function c_2dConnectedComps = get2DConnCompsForEachSlice...
                (x,numCellsR,numCellsC,numSections,neuronIDsForGridCells)

% Extracting 2D connected components (2D neuron sections)
% Give a unique ID for each 2D neuron section

% Output:
%   c_2dConnectedComps - cell array. Each cell containts a cell array as
%   well. Each cell variable of the outer cell array corresponds to an EM slice.
%   Each innter cell variable corresponds to a 2D neuron slice and contains
%   a vector of gridCellIDs 

% Input:
%   neuronIDsForGridCells: N-by-1 matrix where each element says the
%   neuronID the gridCell belongs to.

% Init

