function visualizeNeurons(neuronIDsForGridCells,numCellsR,...
                    numCellsC,numSections,gridSizeR,gridSizeC,rootPixels,...
                    savePath)
                
% assign color according to neuron id

sizeR = gridSizeR * numCellsR;
sizeC = gridSizeC * numCellsC;



% Assign RGB values for each neuronID
uniqueNIDs = unique(neuronIDsForGridCells);
uniqueNIDs = uniqueNIDs(uniqueNIDs>0);
numNeurons = numel(uniqueNIDs);
RGB4NIDs = zeros(numNeurons,3);

for i=1:numNeurons
    R = rand(1); G = rand(1); B = rand(1);
    RGB4NIDs(i,1) = R;
    RGB4NIDs(i,2) = G;
    RGB4NIDs(i,3) = B;    
end

for i=1:numSections
    % visualize each section in RGB
    sectionMat_i = getRGBsection(i,neuronIDsForGridCells,numCellsR,...
                    numCellsC,numSections,gridSizeR,gridSizeC,RGB4NIDs,...
                    sizeR,sizeC,uniqueNIDs,rootPixels);
    figure;imagesc(sectionMat_i);
    % save as png
    filename = sprintf('%d.png',i);
    filename = fullfile(savePath,filename);
    imwrite(sectionMat_i,filename);
end


%% supplementary functions
function sectionMat_i = getRGBsection(sectionID,neuronIDsForGridCells,numCellsR,...
                    numCellsC,numSections,gridSizeR,gridSizeC,RGB4NIDs,...
                    sizeR,sizeC,uniqueNIDs,rootPixels)
%
numGridCellsPerSection = numCellsR * numCellsC;
sectionMat_i = zeros(sizeR,sizeC,3);

Rmat = zeros(sizeR,sizeC);
Gmat = zeros(sizeR,sizeC);
Bmat = zeros(sizeR,sizeC);

gridCellForSection_start = (sectionID-1)*numGridCellsPerSection +1;
gridCellForSection_stop = gridCellForSection_start -1 +...
                            numGridCellsPerSection;
nIDs4gridCells_thisSection = neuronIDsForGridCells...
            (gridCellForSection_start:gridCellForSection_stop);
nIDs_thisSection = unique(nIDs4gridCells_thisSection);
nIDs_thisSection = nIDs_thisSection(nIDs_thisSection>0);
numNIDs_thisSection = numel(nIDs_thisSection);

for i=1:numNIDs_thisSection
    % get the pixels for this neuronID
    nID_i = nIDs_thisSection(i);
    gridCellRelIDs = find(nIDs4gridCells_thisSection==nID_i);
    onCellPixels = getPixIndsForRelGridCellIDsV2...
                        (gridCellRelIDs,rootPixels,gridSizeR,gridSizeC,...
                        sizeR,sizeC);
    % set RGB mat
    nID_logical = (uniqueNIDs==nID_i);
    nID_R = RGB4NIDs(nID_logical,1);
    nID_G = RGB4NIDs(nID_logical,2);
    nID_B = RGB4NIDs(nID_logical,3);
    
    Rmat(onCellPixels) = nID_R;
    Gmat(onCellPixels) = nID_G;
    Bmat(onCellPixels) = nID_B;
    
end

sectionMat_i(:,:,1) = Rmat;
sectionMat_i(:,:,2) = Gmat;
sectionMat_i(:,:,3) = Bmat;