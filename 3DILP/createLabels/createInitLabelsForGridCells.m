function [gridIDs_sectionIDs_rootPixIDsRel,gridCellInteriorlabels,...
            gridCellInteriorRGBLabels,numEdgesY,numEdgesX]...
            = createInitLabelsForGridCells...
            (imageStack3D_label_neuron,numR,numC,numZ,gridResX,gridResY,gridResZ,...
            thresh_mem,imageStack3D_label_mito)

% Input:
%   imageStack3D_label - 3D matrix of the image stack. type: uint8
%       includes membrane cells around the borders
%       
%   thresh_mem = 30; % min percentage of zero pixels to consider the cell to be exterior
% Output:
%   gridIDs_sectionIDs_rootPixIDs - matrix where each col is suggested by name
%   rootPixID is the pixInd of the start pixel (0,0,0) of each cell,
%   wrt each slice (relative coordinates in each slice)
%   gridCellInteriorlabelVector - {0,1}
%       0 - interior
%       1 - exterior (membrane or outer)  
% 

% param 
% thresh_mem = 30; % percentage of zero pixel to consider the cell to be exterior

% grid dimensions (numCells in each dim)

numEdgesX = floor(numC/gridResX);
marginPix_X = mod(numC,gridResX);
gridStartX = floor(marginPix_X/2);
gridStartX = max(1,gridStartX);

numEdgesY = floor(numR/gridResY);
marginPix_Y = mod(numC,gridResY);
gridStartY = floor(marginPix_Y/2);
gridStartY = max(1,gridStartY);

% define node positions
x_pos = gridStartX:gridResX:numC;
y_pos = gridStartY:gridResY:numR;

% Here, a grid node corresponds to the root (0,0,0) of each cube.
% Therefore, we get rid of the final 'grid node' which is the one at the
% boundary
if(x_pos(end)==numC)
    x_pos(end) = [];
end
if(y_pos(end)==numR)
    y_pos(end) = [];
end
[nodeX, nodeY] = meshgrid(x_pos,y_pos);

nodePixList_section = sub2ind([numR numC],nodeY,nodeX);
numGridCellsPerSection = numEdgesX * numEdgesY;
totGridCells = numGridCellsPerSection * numZ;

gridIDs_sectionIDs_rootPixIDsRel = zeros(totGridCells,3);
gridCellInteriorlabels = uint8(zeros(totGridCells,1));
gridCellInteriorRGBLabels = zeros(totGridCells,1);
gridCounter = 0;
for k=1:numZ
    labelImage = imageStack3D_label_neuron(:,:,k);
    % super impose mitochondria on the labels so that the regions
    % indicating mito have RGB val of (0,0,0)
    mitoImage = imageStack3D_label_mito(:,:,k);
    labelImage(mitoImage>0) = 0;
    for j=1:numEdgesX
        for i=1:numEdgesY
            gridCounter = gridCounter + 1;
            % gridCellInd
            gridIDs_sectionIDs_rootPixIDsRel(gridCounter,1) = gridCounter;
            % section_ID
            gridIDs_sectionIDs_rootPixIDsRel(gridCounter,2) = k;
            % rootPixInd
            r0 = nodeY(i,j);
            c0 = nodeX(i,j);
            rootPixInd = nodePixList_section(i,j);
            gridIDs_sectionIDs_rootPixIDsRel(gridCounter,3)...
                = rootPixInd;
            % gridCellLabel
            [label,RGBlabel] = getLabelForCell...
                    (r0,c0,gridResY,gridResX,labelImage,thresh_mem);
            gridCellInteriorlabels(gridCounter) = label;
            gridCellInteriorRGBLabels(gridCounter) = RGBlabel;
            
            
        end
    end
end

%% Supplementary functions
function [label,RGBlabel] = getLabelForCell...
                    (r0,c0,gridResY,gridResX,labelImage_indexed,thresh)
% [sizeR,sizeC] = size(labelImage_indexed);
% get pixels for this cell

% check if it's cellInterior (val>0, label=0) or exterior
% (val=0, label=1). thresh 70% ?
rEnd = r0 + gridResY -1;
cEnd = c0 + gridResX -1;
numPixInCell = gridResY * gridResX;
gridCellPixVal = labelImage_indexed(r0:rEnd,c0:cEnd);
gridCellPixVal = reshape(gridCellPixVal,1,numPixInCell);
% label is 1 (exterior) if either:
%   there are 2 unique nonzero labels or,
%   there are more than 30% zero pixels

% init
label = 0; % default: cellInterior
RGBlabel = mode(gridCellPixVal);
% find non-zero RGB label

uniqueLabels = unique(gridCellPixVal);
numNonZeroLabels = sum(uniqueLabels>0);
if(numNonZeroLabels>1)
    % case where multiple neuron labels are present
    label = 1; % cellExterior
    RGBlabel = 0; 
else
    % percentage of zeropixels
    numZeroPix = sum(sum(gridCellPixVal==0));
    percentageOfZeroPix = numZeroPix*100/(numPixInCell);
    if(percentageOfZeroPix>thresh)
        label = 1; % cellExterior
        RGBlabel = 0;
    end
end





