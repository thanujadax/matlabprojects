function computeFmGridFaces(pathToFm,boundaryGridCellInds,borderFaceInds,...
                gridCIDs_sectionIDs_rootPixIDsRel,numZ,numCellsY,numCellsX,...
                subDir_cellInteriorFm,subDir_cellFaceFm)

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
fm_cellInteriors = readFmFile(pathToGridCellInteriorFeatures,1);
numGridCellFeatures = size(fm_cellInteriors,2);
numCellFaceFeatures = numGridCellFeatures*3;
fm_cellFaces = zeros(numGridCellFacesTot,numCellFaceFeatures);
% numGridCellsPerSlice = numGridCellsTot / numZ;
% cellNeighborMatrix = zeros(numCells,6); % each element is a cellId corresponding to
%                             % the relevant face given by the col number for
%                             % the cellID given by the row number
 
for k=2:numZ-1
    for j=2:numCellsX-1
        for i=2:numCellsY-1
            thisCellInd = sub2ind([numCellsY numCellsX numZ],i,j,k);
            % if this cell is not a boundary cell, get features of all its
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
                fm_cellFaces(faceIndStart:faceIndStop,:)...
                                = FMs_for_directNeighborSet;    
            else 
                disp('computeFmGridFaces.m - Warning 01!')
            
            end % if(sum...) not boundary cell
        end % for i
    end % for j
end % for k
                                                     
% cellFaceIndStop_i = 0;
% for i=1:numGridCellsTot
%     if((sum(boundaryGridCellInds==i))==0)
%         % get ID of direct neighbor for each face in the  given order
% %         cellSectionID = floor(i/numGridCellsPerSlice);
% %         cellID_wrtSection = mod(i,numGridCellsPerSlice);
% %         [cellR,cellC] = ind2sub([numCellsY numCellsX],cellID_wrtSection);
%         [cellR,cellC,cellSectionID] = ind2sub...
%                         ([numCellsY numCellsX numZ],i);
%         setOfDirectNeighbors_6 = getDirectFaceNeighborsInOrder...
%                         (i,cellSectionID,cellR,cellC,...
%                         numCellsY,numCellsX,numZ);
%                     
%         % get rid of zero neighbors?
%         
%         FMs_for_directNeighborSet = readFaceFeatures...
%                         (fm_cellInteriors,i,setOfDirectNeighbors_6);
%         
%         cellFaceIndStart_i = cellFaceIndStop_i +1;
%         cellFaceIndStop_i = cellFaceIndStop_i +6;
%         fm_cellFaces(cellFaceIndStart_i:cellFaceIndStop_i,:)...
%                         = FMs_for_directNeighborSet;
%     else
%         % boundary cell
%         
%     end
%     
% end
%% Saving
% ignore boundary faces?

% Features of different types of faces go in different places
% 1. faces xy {1,2,1,2....}
numXYfaces = numGridCellsTot * 2;
seq = 1:numGridCellsTot;
face1IDs = (seq-1)*6 + 1;
face2IDs = (seq-1)*6 + 2;
face12IDsAll = [face1IDs face2IDs];
% remove border face ids
face12IDsAll = setdiff(face12IDsAll,borderFaceInds);
fm_faces12 = fm_cellFaces(face12IDsAll,:);

fm_name = sprintf('fm_faces12.mat');
saveFileName = fullfile(saveFilePath,fm_name);
save(saveFileName,'fm_faces12');

% 2. faces xz and yz {3,4,5,6,3,4,5,6,...}
seq2 = 1: numGridCellsTot*6;
face3456IDsAll = setdiff(seq2,face12IDsAll);
face3456IDsAll = setdiff(face3456IDsAll,borderFaceInds);
fm_faces3456 = fm_cellFaces(face3456IDsAll,:);

fm_name = sprintf('fm_faces3456.mat');
saveFileName = fullfile(saveFilePath,fm_name);
save(saveFileName,'fm_faces3456');
