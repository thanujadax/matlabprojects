function f = getILP3DGridObjective_traindingLabels...
                    (gridCellInteriorInitLabels,borderCellInds)

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

f(onVarIDs) = reward;
f(offVarIDs) = -reward; % penalty

% borderCell states should be 1 (on) => reward on state
f(borderVarIDs) = borderCellReward;