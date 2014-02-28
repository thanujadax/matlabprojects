function unaryScoresMat = computeUnaryScoreMatRFC...
            (pathToFms,pathToRFCs,numTrees,numCells,NUM_VAR_PER_CELL)
        
% Inputs:
%   RFC - trained random forest classifier

unaryScoresMat = zeros(numCells,NUM_VAR_PER_CELL);
% Each column contains unary scores for the following type of variables
% 1 - cell internal state
% 2 - face xy front
% 3 - face xy back
% 4 - face yz left
% 5 - face yz right
% 6 - face xz top
% 7 - face xz bottom

% load RFC
% load fm
unaryScoresMat(:,1) = getRFCprob(fm,RFC,numTrees); % gridCellProbs_interior
% clear RFC
% clear fm

% load RFC
% load fm
unaryScoresMat(:,2) = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_xy_1 
% clear fm
% load fm
unaryScoresMat(:,3) = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_xy_2 
% clear RFC
% clear fm

% load RFC
% load fm
unaryScoresMat(:,4) = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_yz_3 
% clear fm

% load fm
unaryScoresMat(:,5) = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_yz_4 
% clear fm

% load fm
unaryScoresMat(:,6) = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_xz_5 
% clear fm

% load fm
unaryScoresMat(:,7) = getRFCprob(fm,RFC,numTrees); % gridCellFaceProbs_xz_6 
% clear fm
% clear RFC