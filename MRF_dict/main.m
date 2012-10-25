% MRF learning
% Thanuja

%% parameters, variables and other inputs
hasInitDict = 1;                                % 1 if initial dictionary is provided
pathToDict = 'Dictionary_256x_400w_16bb_NN.mat';
bb = 16;                                        % block size
% get the dictionary words
if(hasInitDict)
    load(pathToDict);       % loads the structure output (from KSVD)
    Dictionary = output.D;  % copy the dictionary
    clear output;           % output from KSVD is no more required 
else
    % generateDictionary()
end
numWords = size(Dictionary,2);

% image for training
trImageFilePath = '/home/thanuja/matlabprojects/data/mitoData/';
imageName = 'stem1_256by256.png';

% Get input image
[IMin,pp]=imread(strcat([trImageFilePath,imageName]));

% imSize = size(IMin_0);
% resizeDim = [imSize(1)-removedEdgeSize imSize(2)-removedEdgeSize];
% IMin = IMin_0(1:resizeDim(1), 1:resizeDim(2), 1);

% fixing input image format
IMin=im2double(IMin);
if (length(size(IMin))>2)
    IMin = rgb2gray(IMin);
end
if (max(IMin(:))<2)
    IMin = IMin*255;
end

% build the association matrices 
% HorizontalAssociations = getHorizontalAssociations(image,Dictionary);
% (remember to add 1)

%% Build MRF
imgInPatches = img2patches(IMin,bb);                        % TODO
% imgInPatches is 3D matrix mxnxk
%   m: number of patches in each column
%   n: number of patches in each row
%   k: number of pixels in each patch = bb*bb

% make edge structure
NN1 = size(imgInPatches,1);
NN2 = size(imgInPatches,2);
numNodes = prod([NN1,NN2]-bb+1);

adj = getAdjMatrix(size(IMin,1),size(IMin,2),bb);           % adjacency matrix
adj = adj + adj';

edgeStruct = UGM_makeEdgeStruct(adj,numWords);
nEdges = edgeStruct.nEdges;

%% Training Ising edge potentials

% Make node map
nodeMap = zeros(numNodes,numWords,'int32');
nodeMap(:,1) = 1;               % initialize with all nodes with word #1. TODO ############

% Make edge map
edgeMap = zeros(numWords,numWords,nEdges,'int32');
% initialize edge map ############## TODO ###################3
for i = 1:numWords
    edgeMap(i,i,:) = numWords;
end

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Example of making potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

% Compute sufficient statistics
suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);

% Evaluate NLL
nll = UGM_MRF_NLL(w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Optimize
w = minFunc(@UGM_MRF_NLL,w,[],nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
nodePot(1,:)
edgePot(:,:,1)
fprintf('(paused)\n');
pause

%% Training (Full Potentials)

edgeMap(2,2,:) = 0;
edgeMap(1,2,:) = 3;
edgeMap(2,1,:) = 4;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Example of making potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

% Compute sufficient statistics
suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);

% Optimize
w = minFunc(@UGM_MRF_NLL,w,[],nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
nodePot(1,:)
edgePot(:,:,1)
fprintf('(paused)\n');
pause


%% Training (incorporating untied potentials for boundaries)

nodeMap(1,1) = 5;
nodeMap(end,1) = 5;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Compute sufficient statistics
suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);

% Optimize
w = minFunc(@UGM_MRF_NLL,w,[],nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
nodePot(1,:)
nodePot(2,:)
edgePot(:,:,1)
fprintf('(paused)\n');
pause

%% Do decoding/infence/sampling in learned model

decode = UGM_Decode_Chain(nodePot,edgePot,edgeStruct)

[nodeBel,edgeBel,logZ] = UGM_Infer_Chain(nodePot,edgePot,edgeStruct);
nodeBel

samples = UGM_Sample_Chain(nodePot,edgePot,edgeStruct);
figure;
imagesc(samples')
title('Samples from MRF model');
fprintf('(paused)\n');
pause

%% Do conditional decoding/inference/sampling in learned model

clamped = zeros(nNodes,1);
clamped(1:2) = 2;

condDecode = UGM_Decode_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Decode_Chain)
condNodeBel = UGM_Infer_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Infer_Chain)
condSamples = UGM_Sample_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Sample_Chain);

figure;
imagesc(condSamples')
title('Conditional samples from MRF model');
fprintf('(paused)\n');
pause



    







