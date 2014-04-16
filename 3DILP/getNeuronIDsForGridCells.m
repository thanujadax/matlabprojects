function [neuronIDsForGridCells, neuronID2DForGridCells] ...
            = getNeuronIDsForGridCells(x,numCellsR,numCellsC,numSections)
    
% Extracting connected components (neurons)

% Outputs:
%   neuronIDsForGridCells - neuronID 0 is for neuron exterior
%   c_gridCellIDsForNeuronIDs

state_neuronInterior = 0;

numGridCells = numCellsR*numCellsC*numSections;

neuronIDsForGridCells = zeros(numGridCells,1);
neuronID2DForGridCells = zeros(numGridCells,1); 

lastUsedID = 0;
freedLabels = [];

lastUsedID_2D = 0;
freedLabels_2D = [];

for k=2:numSections
    for j=2:numCellsC
        for i=2:numCellsR
            
            thisCellInd = sub2ind([numCellsR numCellsC numSections],i,j,k);
            thisCellVarInd = (thisCellInd-1)*7 +1;
            thisCellState = x(thisCellVarInd);
            
            if(thisCellState==state_neuronInterior)
                face1_varInd = thisCellVarInd +1;
                face6_varInd = thisCellVarInd + 6;
                faceStates_thisCell = x(face1_varInd:face6_varInd);

                [setOfDirectNeighbors_6,setOfDirectNeighbors2D_4]...
                        = getDirectFaceNeighborsInOrder...
                                    (thisCellInd,numCellsR,numCellsC,numSections);
                % out of non-membrane neighbors, which are connected across
                % inactive faces
                directNeibhbor_varInds = (setOfDirectNeighbors_6-1)*7 +1;
                directNeibhbor2D_varInds = (setOfDirectNeighbors2D_4-1)*7 +1;
                
                directNeighborStates = x(directNeibhbor_varInds);
                directNeighbor2DStates = x(directNeibhbor2D_varInds);

                neighFaceStates1to6 = getNeighborFaceStates...
                                            (x,directNeibhbor_varInds);
                neighFaceStates3to6 = getNeighborFaceStates...
                                            (x,directNeibhbor2D_varInds);
                nonMemDirectNeighborsCID = getNeighborCIDsToJoin(setOfDirectNeighbors_6,...
                            directNeighborStates,neighFaceStates1to6,...
                            faceStates_thisCell);
                        
                nonMemDirectNeighborsCID_2D = getNeighborCIDsToJoin(setOfDirectNeighbors2D_4,...
                            directNeighbor2DStates,neighFaceStates3to6,...
                            faceStates_thisCell);
                % 3D labels        
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
                        neuronIDsForGridCells(thisCellInd) = nextLabel;
                        [neuronIDsForGridCells,freedLabels] = changeNeighborLabels...
                            (neuronIDsForGridCells,nonMemDirectNeighborsCID,...
                            nextLabel,freedLabels);
                    end

                end
                % 2D labels
                if(numel(nonMemDirectNeighborsCID_2D)>0)
                    
                    nonMemNeighborLabels = neuronID2DForGridCells...
                                    (nonMemDirectNeighborsCID_2D);
                    if(sum(nonMemNeighborLabels)==0)
                        % if they don't have a label, assign a new label to this
                        % gridCell and those neighbors - NEIGHBOR LABELS CHANGED
                        [nextLabel_2D,lastUsedID_2D,freedLabels_2D] ...
                                    = getNextLabel(lastUsedID_2D,freedLabels_2D);
                        neuronID2DForGridCells(thisCellInd) = nextLabel_2D;
                        [neuronID2DForGridCells,freedLabels_2D] = changeNeighborLabels...
                            (neuronID2DForGridCells,nonMemDirectNeighborsCID_2D,...
                            nextLabel_2D,freedLabels_2D);
                        
                    elseif(sum(nonMemNeighborLabels)==1)                        
                        % if one of those neighbors have a label, assign it to all
                        % those neighbors and this cell as well - NEIGHBOR LABELS CHANGED
                        nextLabel_2D = nonMemNeighborLabels(nonMemNeighborLabels>0);
                        neuronID2DForGridCells(thisCellInd) = nextLabel_2D;
                        
                        [neuronID2DForGridCells,freedLabels_2D] = changeNeighborLabels...
                            (neuronID2DForGridCells,nonMemDirectNeighborsCID_2D,...
                            nextLabel_2D,freedLabels_2D);
                        
                    elseif(sum(nonMemNeighborLabels)>1)
                        % if more than one of those neighbors have different
                        % labels, reassign them the same label - NEIGHBOR LABELS CHANGED
                        neighborLabelsAll = nonMemNeighborLabels(nonMemNeighborLabels>0);
                        nextLabel_2D = neighborLabelsAll(1);
                        neuronID2DForGridCells(thisCellInd) = nextLabel_2D;
                        [neuronID2DForGridCells,freedLabels_2D] = changeNeighborLabels...
                            (neuronID2DForGridCells,nonMemDirectNeighborsCID_2D,...
                            nextLabel_2D,freedLabels_2D);
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
        if(label_this~=0)
            cellsWithThisLabel_logicalInds = (neuronIDsForGridCells==label_this);
            neuronIDsForGridCells(cellsWithThisLabel_logicalInds) = newLabel;
            freedLabels = [freedLabels; label_this];
        else
            neuronIDsForGridCells(neighborInds(i)) = newLabel;
            
        end
    end
end

function neighFaceStates1to6 = getNeighborFaceStates(x,neighVarInds)
neighFaceStates1to6 = zeros(6,1);
% for each neighbor in order, get the relevant face ind
% get the states
% for i=1:6
%     neighFaceStates1to6(i) = x(neighVarInds(i)+i);
% end
% face 2 of neighbor 1
neighFaceStates1to6(1) = x(neighVarInds(1)+2);
% face 1 of neighbor 2
neighFaceStates1to6(2) = x(neighVarInds(2)+1);
% face 4 of neighbor 3
neighFaceStates1to6(3) = x(neighVarInds(3)+4);
% face 3 of neighbor 4
neighFaceStates1to6(4) = x(neighVarInds(4)+3);
% face 6 of neighbor 5
neighFaceStates1to6(5) = x(neighVarInds(5)+6);
% face 5 of neighbor 6
neighFaceStates1to6(6) = x(neighVarInds(6)+5);


function neighborsToJoinCIDs = getNeighborCIDsToJoin(neighborCIDs,neighborStates,...
                        neighborFace1to6States,thisCellFaceStates)

zeroStateNeigh_logical = (neighborStates==0);
zeroStateThisCellFaces_logical = (thisCellFaceStates==0);
zeroStateNeighborFaces1to6_logical = (neighborFace1to6States==0);

neighborsToJoin_logical = zeroStateNeigh_logical & ...
                    zeroStateThisCellFaces_logical & ...
                    zeroStateNeighborFaces1to6_logical;
                
neighborsToJoinCIDs = neighborCIDs(neighborsToJoin_logical);
