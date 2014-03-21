function neuronIDsForGridCells ...
            = getNeuronIDsForGridCells(x,numCellsR,numCellsC,numSections)
    
% Extracting connected components (neurons)

% Outputs:
%   neuronIDsForGridCells - neuronID 0 is for neuron exterior
%   c_gridCellIDsForNeuronIDs

state_neuronInterior = 0;

numGridCells = numCellsR*numCellsC*numSections;

neuronIDsForGridCells = zeros(numGridCells,1);
lastUsedID = 0;
freedLabels = [];

for k=2:numCells
    for j=2:numCellsC
        for i=2:numCellsR
            
            flag_neighborsChanged = 0;
            
            thisCellInd = sub2ind([numCellsR numCellsC numSections],i,j,k);
            thisCellVarInd = (thisCellInd-1)*7 +1;
            thisCellState = x(thisCellVarInd);
            
            if(thisCellState==state_neuronInterior)
                face1_varInd = thisCellVarInd +1;
                face6_varInd = thisCellVarInd + 6;
                faceStates_thisCell = x(face1_varInd:face6_varInd);

                setOfDirectNeighbors_6 = getDirectFaceNeighborsInOrder...
                                    (thisCellInd,numCellsR,numCellsC,numSections);
                % out of non-membrane neighbors, which are connected across
                % inactive faces
                directNeibhbor_varInds = (setOfDirectNeighbors_6-1)*7 +1; 
                directNeighborStates = x(directNeibhbor_varInds);
                nonMemDirectNeighborsCID = setOfDirectNeighbors_6(...
                            directNeighborStates==state_neuronInterior);
                        
                if(numel(nonMemDirectNeighborsCID)>0)
                    
                    nonMemNeighborLabels = neuronIDsForGridCells...
                                    (nonMemDirectNeighborsCID);
                                
                    if(sum(nonMemNeighborLabels)==0)
                        % if they don't have a label, assign a new label to this
                        % gridCell and those neighbors - NEIGHBOR LABELS CHANGED
                        [nextLabel,lastUsedID,freedLabels] ...
                                    = getNextLabel(lastUsedID,freedLabels);
                        neuronIDsForGridCells(thisCellInd) = nextLabel;
                        [neuronIDsForGridCells,freedLabels] = changeNeighborLabels...
                            (neuronIDsForGridCells,nonMemDirectNeighborsCID,...
                            nextLabel,freedLabels);
                        
                    elseif(sum(nonMemNeighborLabels)==1)                        
                        % if one of those neighbors have a label, assign it to all
                        % those neighbors and this cell as well - NEIGHBOR LABELS CHANGED
                        nextLabel = nonMemNeighborLabels(nonMemNeighborLabels>0);
                        neuronIDsForGridCells(thisCellInd) = nextLabel;
                        
                        [neuronIDsForGridCells,freedLabels] = changeNeighborLabels...
                            (neuronIDsForGridCells,nonMemDirectNeighborsCID,...
                            nextLabel,freedLabels);
                        
                    elseif(sum(nonMemNeighborLabels)>1)
                        % if more than one of those neighbors have different
                        % labels, reassign them the same label - NEIGHBOR LABELS CHANGED
                        neighborLabelsAll = nonMemNeighborLabels(nonMemNeighborLabels>0);
                        nextLabel = neighborLabelsAll(1);
                        [neuronIDsForGridCells,freedLabels] = changeNeighborLabels...
                            (neuronIDsForGridCells,nonMemDirectNeighborsCID,...
                            nextLabel,freedLabels);
                    end

                end
            end
        end
    end
end



%% Supplementary functions
function [nextLabel,lastUsedID,freedLabels] = getNextLabel(lastUsedID,freedLabels)
    if(numel(freedLabels)>0)
        nextLabel = freedLabels(1);
        freedLabels(1) = [];
    else
        nextLabel = lastUsedID + 1;
        lastUsedID = nextLabel;
        
    end
    
function [neuronIDsForGridCells,freedLabels] = changeNeighborLabels...
                        (neuronIDsForGridCells,neighborInds,newLabel,freedLabels)
% 
numNeighbors = numel(neighborInds);

for i=1:numNeighbors
    % get the label for this neighbor
    label_this = neuronIDsForGridCells(neighborInds(i));
    if(label_this~=newLabel)
        % if it's different from the newLabel, get all the cells with this
        % label and assign them the newLabel.
        cellsWithThisLabel_logicalInds = (neuronIDsForGridCells==label_this);
        neuronIDsForGridCells(cellsWithThisLabel_logicalInds) = newLabel;
        freedLabels = [freedLabels; label_this];
    end
end
