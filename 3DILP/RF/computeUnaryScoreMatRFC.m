function unaryScoresMat = computeUnaryScoreMatRFC...
            (pathToFms,pathToRFCs,numTrees,numCells,NUM_VAR_PER_CELL,...
            borderCellFaceInds)
        
% Inputs:
%   RFC - trained random forest classifier

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

fm_gridCells = 'fm_cellInterior.mat'; % all cells including border cells
fm_faces12 = 'fm_faces12.mat';
fm_faces34 = 'fm_faces34.mat';
fm_faces56 = 'fm_faces56.mat';

%% gridCellProbs_interior
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_gridCells);
RFC = importdata(RFC_filename);
% load fm
fm_filename = fullfile(pathToFms,fm_gridCells);
fm = importdata(fm_filename);
cellInteriorScores = getRFCprob(fm,RFC,numTrees); 
clear RFC fm
%% faces 1,2
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_faces12);
RFC = importdata(RFC_filename);
% load fm
fm_filename = fullfile(pathToFms,fm_faces12);
fm = importdata(fm_filename);
face12Scores = getRFCprob(fm,RFC,numTrees);
clear RFC fm
%% faces 3,4
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_faces3456);
RFC = importdata(RFC_filename);
% load fm
fm_filename = fullfile(pathToFms,fm_faces34);
fm = importdata(fm_filename);
face34Scores = getRFCprob(fm,RFC,numTrees); 
clear RFC fm 
%% faces 5,6
% load RFC
RFC_filename = fullfile(pathToRFCs,forest_faces3456);
RFC = importdata(RFC_filename);
% load fm
fm_filename = fullfile(pathToFms,fm_faces56);
fm = importdata(fm_filename);
face56Scores = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_xy_1 
clear RFC fm
%% Create unary score matrix for all ILP variables
% {cell1,face1,face2,f3,f4,f5,f6} 
% {cell2...} ...
if(NUM_VAR_PER_CELL==6)
    unaryScoresMat = zeros(numCells,NUM_VAR_PER_CELL);
    % col1: gridCellScores
    unaryScoresMat(:,1) = cellInteriorScores;
    
    % col2: face1
    face1LInd = 1:2:numCells*2;
    face1LInd = setdiff(face1LInd,borderCellFaceInds);
    face1Probs = zeros(numCells,1);
    face12_seq1 = 1:2:numel(face12Scores);
    face1Probs(face1LInd) = face12Scores(face12_seq1);
    unaryScoresMat(:,2) = face1Probs;
    
    % col3: face2
    face2LInd = 2:2:numCells*2;
    face2LInd = setdiff(face2LInd,borderCellFaceInds);
    face2Probs = zeros(numCells,1);
    face12_seq2 = 2:2:numel(face12Scores);
    face2Probs(face2LInd) = face12Scores(face12_seq1);
    unaryScoresMat(:,3) = face2Probs;
    
else
    unaryScoreMat = 0;
    disp('Error: computeUnaryScoreMatRFC')
end
