function [mClusters,nCC] = createAllClusters(cImages,varargin)
% Create one cluster for all the membranes (all the orientations & 
% junctions), one cluster per neuron id (contaning the mitochondria), 
% one cluster per glia id, one cluster per synapse id.

% In our setting
if nargin < 3
    vMembranesMergeId = [1:1:5];
end
if nargin < 4
    vNeuronsMergeId = [7,9];
end
if nargin < 5
    vGliaMergeId = [6];
end
if nargin < 5
    vSynapseMergeId = [8];
end

% Initializations
mClusters = zeros(size(cImages{1}));
nNCC = 4;
nCC= 0;

% One cluster for the membranes (all the orientations & junctions)
for i=vMembranesMergeId
    mClusters(cImages{i} ~= 0) = 1;
end
nCC = nCC + 1;

% One cluster per neuron id (containing the mitochondria)
mTemp = zeros(size(cImages{1}));
for i=vNeuronsMergeId
    mTemp(cImages{i} ~= 0) = 1;
end
[mCC,nCCTemp] = bwlabel(mTemp,nNCC);
mClusters(mTemp ~= 0) = mCC(mTemp ~= 0) + nCC;

nCC = nCC + nCCTemp;

% One cluster per glia id
mTemp = zeros(size(cImages{1}));
for i=vGliaMergeId
    mTemp(cImages{i} ~= 0) = 1;
end
[mCC,nCCTemp] = bwlabel(mTemp,nNCC);

mClusters(mTemp ~= 0) = mCC(mTemp ~= 0) + nCC;
nCC = nCC + nCCTemp;

% One cluster per synapse id
mTemp = zeros(size(cImages{1}));
for i=vSynapseMergeId
    mTemp(cImages{i} ~= 0) = 1;
end
[mCC,nCCTemp] = bwlabel(mTemp,nNCC);

mClusters(mTemp ~= 0) = mCC(mTemp ~= 0) + nCC;
nCC = nCC + nCCTemp;

function setdefaults(defaults)
    var_names = evalin('caller', 'whos');
    default_names = fieldnames(defaults);
    must_replace = ~ismember(default_names, var_names);
    for k = find(must_replace)
        name = default_names(k);
        assignin('caller', name, defaults.(name));
    end

