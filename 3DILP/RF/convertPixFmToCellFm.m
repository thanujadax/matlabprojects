function fm_gridsOfSlice_i = convertPixFmToCellFm...
        (fm_slice,rootPixIndsOfCells,gridResX,gridResY,sizeR,sizeC,...
        numGridCellsPerSection,numGridCellFeatures)

% Using the pixel feature matrix of one slice, compute the feature matrix
% for gridCells

% Inputs:
%   fm_slice - featureMatrix for the current slice (3D)
%   gridCIDs_

% The input section corresponds to a set of gridCells

% For each gridCell, get the set of pixInds
% -> get the relevant entries from fm_slice

[rStart,cStart] = ind2sub([sizeR sizeC],rootPixIndsOfCells);

rStop = rStart + gridResY -1;
cStop = cStart + gridResX -1;

% get pixels
% numPixPerCell = gridResX * gridResY;

% Init
fm_gridsOfSlice_i = zeros(numGridCellsPerSection,numGridCellFeatures);


% compute relevant grid features (mean,min,max,var,median)
for i=1:numGridCellsPerSection
    fm_slice_tmp = fm_slice;
    rCell = rStart(i) : rStop(i);
    cCell = cStart(i) : cStop(i); 
    
%     pixIndsCell = sub2ind([sizeR sizeC],rCell,cCell);
    
    fm_pixCell = fm_slice_tmp(rCell,cCell,:);
    
    fm_gridsOfSlice_i(i,:) = computeGridFeaturesFromPixFeatures...
                (fm_pixCell);
    
end