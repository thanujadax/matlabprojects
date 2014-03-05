function section = visualizationOnGridCells_singleSection(gridCellStates,rootPixels,...
                        gridResY,gridResX,sizeR,sizeC)

% Input:
%   gridCellStates - vector containing the states of the gridCells in one
%   section.

% mapping from gridCellID(relative) to set of pixels in the section

% gridCell 

section = zeros(sizeR,sizeC);

% get active grids so that the 
interiorCellPixels = getPixIndsForRelGridCellIDs...
                        (gridCellStates,rootPixels,gridResY,gridResX,...
                        sizeR,sizeC);
section(interiorCellPixels) = 1;
