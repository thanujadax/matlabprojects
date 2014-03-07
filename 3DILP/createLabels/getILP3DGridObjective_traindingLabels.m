function f = getILP3DGridObjective_traindingLabels...
                    (gridCellInteriorInitLabels,borderCellInds,...
                    gridCellFaceInitLabels)

% assign a reward for agreeing with the input initial labels

% param
reward = -5;
borderCellReward = -1000;

numCells = numel(gridCellInteriorInitLabels);
totNumVar = numCells * 7;

% init
f = zeros(totNumVar,1);

onCellIDs = find(gridCellInteriorInitLabels>0);
offCellIDs = find(gridCellInteriorInitLabels==0);

onVarIDs = (onCellIDs-1) .*7 +1;
offVarIDs = (offCellIDs-1) .*7 +1;
borderVarIDs = (borderCellInds-1) .*7 +1;

% f(onVarIDs) = reward;
% f(offVarIDs) = -reward; % penalty


% Create matrix where each row corresponds to one gridCell. 1st col is the
% gridCellInterior state. cols 2:7 correspond to the states of the faces

stateMatrix = zeros(numCells,7);

stateMatrix(:,1) = gridCellInteriorInitLabels;
stateMatrix(:,2:7) = gridCellFaceInitLabels;

stateMatrix = stateMatrix'; % now each col corresponds to one gridCell

initVariableStates = repmat(stateMatrix,1,totNumVar);

f(initVariableStates) = reward;
f(~initVariableStates) = -reward;

% borderCell states should be 1 (on) => reward on state
f(borderVarIDs) = borderCellReward;