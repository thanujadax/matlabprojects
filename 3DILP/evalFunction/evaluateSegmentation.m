% evaluation of 3D segmentation
function s_evaluation = evaluateSegmentation()
% Inputs
%   Initial 3D segmentation
%   Trained classifiers
%   Data ? ()

% Outputs
%   Global constraint violations
%       mutually exclusive sets of links
%   Abrupt ends
%   locally confident structures (medium range sub-graphs)
%   weak medium range sub-graphs
%   weak 3D links
%   weak (2D) nodes
%       potential merge errors
%       potential split errors


%% Extract components, and their features
% 2DnodeIDs = {sectionID,rID,cID,centrOfMass(pixInd_inSection),neuronID}
% 3DshortRangeComponents = {numFacesSect1,numfacesSect2,[2DnodeIDs]}


[neuronIDsForGridCells, nodeIDs2DForGridCells] ...
            = getNeuronIDsForGridCells(x,numCellsR,numCellsC,numSections);

neuronComp2Dfeatures = getNeuron2Dfeatures();
neuronComp3Dfeatures = getNeuron3Dfeatures();
skeletonFeatures = getSkeletonFeatures();

%%  Evaluate components

