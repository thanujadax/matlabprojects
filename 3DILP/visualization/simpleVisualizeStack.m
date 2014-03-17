function simpleVisualizeStack(x,rootPixels,gridResY,gridResX,sizeR,sizeC,...
                    numCellsR,numCellsC,numSections)

% Inputs
%   x
%   rootPixels: pixInds of gridCells in order wrt a single section


% Shows all cell interior as bright (1) and exterior as dark (0)
% for all the sections in the input stack

% separate indicator variables for each stack
% separate the indicator variables for cellAcitivation
% starting from the 1st, every 7th variable is a gridCellState variable

totVar = numel(x);
cellStateListInds = 1:7:totVar;
gridCellStateVariables = x(cellStateListInds);

numCellsPerSection = numCellsR * numCellsC;

gridCellStateColumns = reshape(gridCellStateVariables,...
                    numCellsPerSection,numSections);

% call visualize_section_simple for each section


for i=1:numSections
    gridCellStates = gridCellStateColumns(:,i);
    section = visualizationOnGridCells_singleSection(gridCellStates,rootPixels,...
                        gridResY,gridResX,sizeR,sizeC);
    title_str = sprintf('ILP segmentation cell interior: section %d',i);
    figure;imshow(section);title(title_str)
end