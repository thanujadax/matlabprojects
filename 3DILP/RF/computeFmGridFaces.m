function computeFmGridFaces(pathToFm,boundaryGridCellInds,borderFaceInds,...
                gridCIDs_sectionIDs_rootPixIDsRel,numZ,numCellsY,numCellsX,...
                subDir_cellInteriorFm,subDir_cellFaceFm,...
    listInds_fm_face12_name,listInds_fm_face34_name,listInds_fm_face56_name,...
    name_fm_faces12,name_fm_faces34,name_fm_faces56,name_fm_cellInterior)

% Inputs:

% pathToFm - directory in which the subdirectories containing fm files are
% located
% subDir_cellInteriorFm - name of subdirectory in which the feature mats of
% cell interiors are already saved (inside pathToFm)
% subDir_cellFaceFm - name of subdirectory in which the feature mats of
% cell faces are to be saved (inside pathToFm)

% Outputs:
% no value is returned. The following two files are saved in the specified
% output path (pathToFm/subDir_cellFaceFm)
%   fm_faces12.mat
%   fm_faces3456.mat          

% computes features for each face of eah grid cell
% takes into account the features of the direct neighbor of the face as
% well as the other neighbors around it.

% FEATURES ARE NOT CALCULATED FOR GRID BOUNDARY FACES!

% version 1.0

% TODO: version 1.0 takes into account only the direct neighbor of the
% face. Extend it to consider the features of other neighbors as well.

% each grid cell has 6 faces
% should not consider the boundary cells

% Face order
%   xy1(front),xy2(back),yz1(left),yz2(right),xz1(top),xz2(bottom)

% for each cell
% get the neighbors in order (starting with the direct neighbor)
% calculate differential features (direct neighbor)
% combine all features to one vector 
% save all faces in one file

% __________________
% |       |         |
% |       |         |
% |  A   a|b    B   |
% |       |         |
% |_______|_________|

% feature vector for face a = {f(A-B),f(A),f(B)}

% path to save section features
saveFilePath = fullfile(pathToFm,subDir_cellFaceFm);

% path to read already calculated cell interior feature mat
pathToGridCellInteriorFeatures = fullfile(pathToFm,subDir_cellInteriorFm);

% check if subdirectory exists. Create if not
checkAndCreateSubDir(pathToFm,subDir_cellFaceFm);

% featureFile.mat name structure: fm_slice_%d.mat

numGridCellsTot = size(gridCIDs_sectionIDs_rootPixIDsRel,1); % including boundary cells
numGridCellFacesTot = numGridCellsTot * 6;

% Init
% fm_cellInteriors = readFmFile(pathToGridCellInteriorFeatures,1);
fm_cellInteriorFileName = fullfile(pathToGridCellInteriorFeatures,name_fm_cellInterior);
fm_cellInteriors = importdata(fm_cellInteriorFileName); % with border cells
numGridCellFeatures = size(fm_cellInteriors,2);
numCellFaceFeatures = numGridCellFeatures*3;
fm_cellFaces = zeros(numGridCellFacesTot,numCellFaceFeatures);
% numGridCellsPerSlice = numGridCellsTot / numZ;
% cellNeighborMatrix = zeros(numCells,6); % each element is a cellId corresponding to
%                             % the relevant face given by the col number for
%                             % the cellID given by the row number
numNonBoundaryCells = (numCellsY-2)*(numCellsX-2)*(numZ-2);
numNonBoundaryCellFaces = numNonBoundaryCells * 6;
faceIndsList_all = zeros(numNonBoundaryCellFaces,1); 
listStop = 0;
for k=2:numZ-1
    for j=2:numCellsX-1
        for i=2:numCellsY-1
            thisCellInd = sub2ind([numCellsY numCellsX numZ],i,j,k);
            % this cell is not a boundary cell. Get features of all its
            % six faces. features of borderCellFaces are automatically kept
            % zero due to init and are discarded when saving feature mat.
            if(sum(intersect(boundaryGridCellInds,thisCellInd))==0)
                % get the neighbor cells
                setOfDirectNeighbors_6 = getDirectFaceNeighborsInOrder...
                                (thisCellInd,numCellsY,numCellsX,numZ);                           
                FMs_for_directNeighborSet = readFaceFeatures...
                        (fm_cellInteriors,thisCellInd,setOfDirectNeighbors_6);
                faceIndStart = (thisCellInd-1)*6 +1;
                faceIndStop = faceIndStart +5;
                listStart = listStop +1;
                listStop = listStop +6;
                faceIndsList_all(listStart:listStop)...
                    =  faceIndStart:faceIndStop;
                fm_cellFaces(faceIndStart:faceIndStop,:)...
                                = FMs_for_directNeighborSet;    
            else 
                disp('computeFmGridFaces.m - Warning 01!')
            
            end % if(sum...) not boundary cell
        end % for i
    end % for j
end % for k
                                                     
%% Saving
% ignore boundary faces?

% Features of different types of faces go in different places
% 1. faces xy {1,2,1,2....}
face1LIDs = 1:6:numNonBoundaryCellFaces;
face2LIDs = 2:6:numNonBoundaryCellFaces;
listInds_fm_face12 = [face1LIDs face2LIDs];
listInds_fm_face12 = sort(listInds_fm_face12);
listInds_fm_face12 = faceIndsList_all(listInds_fm_face12,:);
fm_faces12 = fm_cellFaces(listInds_fm_face12,:);
% save face12 fm
% fm_name = sprintf('fm_faces12.mat');
saveFileName = fullfile(saveFilePath,name_fm_faces12);
save(saveFileName,'fm_faces12');
% save face12 faceInds
saveFileName = fullfile(saveFilePath,listInds_fm_face12_name);
save(saveFileName,'listInds_fm_face12');

clear face1IDs face2IDs fm_faces12 face12IDsAll

% 2. faces xz {3,4,3,4,...}
face3LIDs = 3:6:numNonBoundaryCellFaces;
face4LIDs = 4:6:numNonBoundaryCellFaces;
listInds_fm_face34 = [face3LIDs face4LIDs];
listInds_fm_face34 = sort(listInds_fm_face34);
listInds_fm_face34 = faceIndsList_all(listInds_fm_face34,:);
fm_faces34 = fm_cellFaces(listInds_fm_face34,:);
% save face34 fm
% fm_name = sprintf('fm_faces34.mat');
saveFileName = fullfile(saveFilePath,name_fm_faces34);
save(saveFileName,'fm_faces34');
% save face34 faceInds
saveFileName = fullfile(saveFilePath,listInds_fm_face34_name);
save(saveFileName,'listInds_fm_face34');

clear face3IDs face4IDs fm_faces34 face34IDsAll

% 3. faces yz {5,6,5,6,...}
face5LIDs = 5:6:numNonBoundaryCellFaces;
face6LIDs = 6:6:numNonBoundaryCellFaces;
listInds_fm_face56 = [face5LIDs face6LIDs];
listInds_fm_face56 = sort(listInds_fm_face56);
listInds_fm_face56 = faceIndsList_all(listInds_fm_face56,:);
fm_faces56 = fm_cellFaces(listInds_fm_face56,:);
% save face56 fm
% fm_name = sprintf('fm_faces56.mat');
saveFileName = fullfile(saveFilePath,name_fm_faces56);
save(saveFileName,'fm_faces56');
% save face56 faceInds
saveFileName = fullfile(saveFilePath,listInds_fm_face56_name);
save(saveFileName,'listInds_fm_face56');
clear face5IDs face6IDs fm_faces56 face56IDsAll

