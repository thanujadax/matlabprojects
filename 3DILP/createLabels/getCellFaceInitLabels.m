function gridCellFaceInitLabels = getCellFaceInitLabels...
                    (cellInteriorRGBLabelsAll,numSections,...
                    numEdgesX,numEdgesY)

% Inputs:
%   cellInteriorRGBLabelsAll - unique RGB labels of annotated neurons
%       0 for cellExterior
%       RGB value (int) for cell interior (majority vote for larger )
%   nodeX (=nodeC) -  all the x components of pixelInds for gridCell
%   rootPixels
%   nodeY (=nodeR) -  all the y components of pixelInds for gridCell
%   rootPixels

% For each cell take its 6 faces (neighbors) in order and compare the
% labels thisCell and neighborCell to get the label for the corresponding
% face.

% only have to check faces 2,4 and 6. By symmetry, the faces 1,3,5 of the
% neighbors can be labeled as well, after the same comparison.


numFacesPerCell = 6;
totNumCells = numel(cellInteriorRGBLabelsAll);
gridCellFaceInitLabels = zeros(totNumCells,numFacesPerCell);

for k=1:numSections
    for j=1:numEdgesX
        for i=1:numEdgesY
            % debug
            if(i==2 && j==2 && k==2)
                aa=99;
            end
            %
           
            thisCellInd = sub2ind([numEdgesY numEdgesX numSections],i,j,k);
            setOfDirectNeighbors_3 = getDirectFaceNeighborsInOrder_246...
                                (thisCellInd,k,i,j,...
                                numEdgesY,numEdgesX,numSections);
            % contains zero when there's no neighbor!!!!!!!!!!
            setOfDirectNeighbors_3 = setOfDirectNeighbors_3...
                                    (setOfDirectNeighbors_3>0);
            
            thisCellLabel = cellInteriorRGBLabelsAll(thisCellInd);
            neighbor246Labels = cellInteriorRGBLabelsAll(setOfDirectNeighbors_3);
            
            numNeighbors = numel(setOfDirectNeighbors_3);          
            faceLabels246_135 = getFaceLabels246_135...
                                    (thisCellLabel,neighbor246Labels,numNeighbors);
                        
            if(numNeighbors==3)
                % get labels for faces 2,4,6 of this cell
                % these are the same labels for the faces 1,3,5 of the 2,4,6
                % neighbors respectively.
                
                thisLabelFaceInds = [2 4 6];
                gridCellFaceInitLabels(thisCellInd,thisLabelFaceInds) = ...
                            faceLabels246_135(:,1);
                % neighbor1
                neigh = 1;
                neigh_face = 1;
                gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                            = faceLabels246_135(neigh,2);

                % neighbor2
                neigh = 2;
                neigh_face = 3;
                gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                            = faceLabels246_135(neigh,2);

                % neighbor3
                neigh = 3;
                neigh_face = 5;
                gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                            = faceLabels246_135(neigh,2);
            elseif(numNeighbors==2)
                if(k==numSections)
                    thisLabelFaceInds = [4 6];
                    
                    gridCellFaceInitLabels(thisCellInd,thisLabelFaceInds) = ...
                                faceLabels246_135(:,1);
                    % neighbor1
                    neigh = 1;
                    neigh_face = 3;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);

                    % neighbor2
                    neigh = 2;
                    neigh_face = 5;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);
                elseif(j==numEdgesX)
                    thisLabelFaceInds = [2 6];
                    gridCellFaceInitLabels(thisCellInd,thisLabelFaceInds) = ...
                                faceLabels246_135(:,1);
                    % neighbor1
                    neigh = 1;
                    neigh_face = 1;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);

                    % neighbor2
                    neigh = 2;
                    neigh_face = 5;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);
                    
                elseif(i==numEdgesY)
                    thisLabelFaceInds = [2 4];
                    gridCellFaceInitLabels(thisCellInd,thisLabelFaceInds) = ...
                                faceLabels246_135(:,1);
                    % neighbor1
                    neigh = 1;
                    neigh_face = 1;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);

                    % neighbor2
                    neigh = 2;
                    neigh_face = 3;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);
                end
                
            elseif(numNeighbors==1)
                if(k~=numSections)
                    thisLabelFaceInds = 2;
                    gridCellFaceInitLabels(thisCellInd,thisLabelFaceInds) = ...
                                faceLabels246_135(:,1);
                    % neighbor1
                    neigh = 1;
                    neigh_face = 1;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);
                elseif(j~=numEdgesX)
                    thisLabelFaceInds = 4;
                    gridCellFaceInitLabels(thisCellInd,thisLabelFaceInds) = ...
                                faceLabels246_135(:,1);
                    % neighbor1
                    neigh = 1;
                    neigh_face = 3;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);
                else
                    thisLabelFaceInds = 6;
                    gridCellFaceInitLabels(thisCellInd,thisLabelFaceInds) = ...
                                faceLabels246_135(:,1);
                    % neighbor1
                    neigh = 1;
                    neigh_face = 5;
                    gridCellFaceInitLabels(setOfDirectNeighbors_3(neigh),neigh_face)...
                                = faceLabels246_135(neigh,2);
                end
                
                
            end
            
        end
    end
end

%% Supplementary functions

function faceLabels246_135 = getFaceLabels246_135...
                    (thisCellRGBLabel,neighbor246RGBLabels,numNeighbors)
% Inputs:
%   thisCellLabel
%   neighbor246Labels

% Output:
%   faceLabels246_135 - 3x2 matrix
%       - col1 has the faceStates(2,4,6) of thisCell
%       - col2 has the faceStates(1,3,5) of the corresponding neighbors

% logic:
% thisCellState thisFace thatFace thatCell
%   0           0           0       0
%   0           1           1       0
%   0           1           0       1
%   1           0           0       1

% RGB labels
%   x           0           0       x
%   x           1           1       y
%   x           1           0       0
%   0           0           0       0

%numNeighbors = sum(neighbor246Labels>0);

faceLabels246_135 = zeros(numNeighbors,2);

thisCellLabels_3 = ones(1,numNeighbors) .* thisCellRGBLabel;

tmp1 = neighbor246RGBLabels' - thisCellLabels_3;
% this case is already taken care of by initialization
%   faceLabels246_135(tmp==0,:) = 0;

nzdiff_logicalInd = (tmp1~=0);
if(sum(nzdiff_logicalInd)>0)
    tmp2 = nzdiff_logicalInd & logical(thisCellLabels_3);
    if(sum(tmp2)>0)
        faceLabels246_135(tmp2,1) = 1;
    end
    tmp3 = nzdiff_logicalInd' & logical(neighbor246RGBLabels);
    if(sum(tmp3)>0)
        faceLabels246_135(tmp3,2) = 1;
    end
end
