function unaryScoresMat = computeUnaryScoreMatRFC...
            (pathToFms,pathToRFCs,numTrees,numCells,NUM_VAR_PER_CELL,...
            subDir_cellInteriorFm,subDir_cellFaceFm,...
        name_fm_faces12,listInds_fm_face12_name,...
        name_fm_faces34,listInds_fm_face34_name,...
        name_fm_faces56,listInds_fm_face56_name,...
        name_fm_cellInterior,listInds_fm_cells_name)
        
% Inputs:
%   RFC - trained random forest classifier

% If NUM_VAR_PER_CELL == 7:
% Each column contains unary scores for the following type of variables
% 1 - cell internal state
% 2 - face xy front
% 3 - face xy back
% 4 - face yz left
% 5 - face yz right
% 6 - face xz top
% 7 - face xz bottom

%% File names
forest_gridCells = 'RFC_gridCells.mat';
forest_faces12 = 'RFC_gridFaces12.mat';
forest_faces3456 = 'RFC_gridFaces3456.mat';

% name_fm_cellInterior = 'fm_cellInterior.mat'; % all cells including border cells
% name_fm_faces12 = 'fm_faces12.mat';
% name_fm_faces34 = 'fm_faces34.mat';
% name_fm_faces56 = 'fm_faces56.mat';

%% gridCellProbs_interior
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_gridCells);
RFC = importdata(RFC_filename);
% load fm
fm_CellInteriorDir = fullfile(pathToFms,subDir_cellInteriorFm);
fm_filename = fullfile(fm_CellInteriorDir,name_fm_cellInterior);
fm = importdata(fm_filename);
cellInteriorScores = getRFCprob(fm,RFC,numTrees); 
clear RFC fm
%% faces 1,2
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_faces12);
RFC = importdata(RFC_filename);
% load fm
fm_facesDir = fullfile(pathToFms,subDir_cellFaceFm);
fm_filename = fullfile(fm_facesDir,name_fm_faces12);
fm = importdata(fm_filename);
face12Scores = getRFCprob(fm,RFC,numTrees);
clear RFC fm
%% faces 3,4
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_faces3456);
RFC = importdata(RFC_filename);
% load fm
fm_filename = fullfile(fm_facesDir,name_fm_faces34);
fm = importdata(fm_filename);
face34Scores = getRFCprob(fm,RFC,numTrees); 
clear RFC fm 
%% faces 5,6
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_faces3456);
RFC = importdata(RFC_filename);
% load fm
fm_filename = fullfile(fm_facesDir,name_fm_faces56);
fm = importdata(fm_filename);
face56Scores = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_xy_1 
clear RFC fm
%% Create unary score matrix for all ILP variables
% {cell1,face1,face2,f3,f4,f5,f6} 
% {cell2...} ...
if(NUM_VAR_PER_CELL==7)
    unaryScoresMat = zeros(numCells,NUM_VAR_PER_CELL);
    % col1: gridCellScores
    numCellScores = numel(cellInteriorScores);
    % load cellIndsList_forFM    
    fm_dir_cellInterior = fullfile(pathToFms,subDir_cellInteriorFm);
    cellIndsListFile = fullfile(fm_dir_cellInterior,listInds_fm_cells_name);
    cellIndsList_forFM = importdata(cellIndsListFile);
    cellScoreVector = zeros(numCells,1);
    cellScoreVector(cellIndsList_forFM) = cellInteriorScores;
    
    % cell face matrix (transposed) - filled with face scores for each cell
    faceScoreMatrix_t = zeros(6,numCells);
    % faces 1,2
    fm_dir_faces = fullfile(pathToFms,subDir_cellFaceFm);
    faces12IndsListFile = fullfile(fm_dir_faces,listInds_fm_face12_name);
    faces12Inds_forFM = importdata(faces12IndsListFile);
    faceScoreMatrix_t(faces12Inds_forFM) = face12Scores;
    % faces 3,4
    faces34IndsListFile = fullfile(fm_dir_faces,listInds_fm_face34_name);
    faces34Inds_forFM = importdata(faces34IndsListFile);
    faceScoreMatrix_t(faces34Inds_forFM) = face34Scores;
    % faces 5,6
    faces56IndsListFile = fullfile(fm_dir_faces,listInds_fm_face56_name);
    faces56Inds_forFM = importdata(faces56IndsListFile);
    faceScoreMatrix_t(faces56Inds_forFM) = face56Scores;
    
    % fill in unaryScoreMat
    unaryScoresMat(:,1) = cellScoreVector;
    unaryScoresMat(:,2:7) = faceScoreMatrix_t';
    
else
    unaryScoresMat = 0;
    disp('Error: computeUnaryScoreMatRFC')
end
