function f2 = getPairwisePotentials(patchID,currentLabels,verticalAMat,horizontalAMat,rowSize)
% for a single patch considering its neighborhoot
    pairwisePotMat = zeros(1,4);

    if(patchID - rowSize > 0)
        pairwisePotMat(1,1) = -1.*log(verticalAMat(currentLabels(patchID),currentLables(patchID-rowSize)));
    else
        pairwisePotMat(1,1) = 0;    % nothing on top
    end
    % condCostGivenBottom
    if(patchID+rowSize <= size(currentLabels,2))
        pairwisePotMat(1,2) = -1.*log(verticalAMat(currentLables(patchID+rowSize),currentLabels(patchID)));
    else
        pairwisePotMat(1,2) = 0;    % nothing to the bottom
    end
    % condCostGivenLeft
    if(patchID-1 > 0 && mod(patchID,rowSize)>1)
        pairwisePotMat(1,3) = -1.*log(horizontalAMat(currentLabels(patchID),currentLabels(patchID-1)));
    else
        pairwisePotMat(1,3) = 0;    % nothing to the left
    end
    % condCostGivenRight
    if(patchID+1<=size(currentLabels,2) && mod(patchID,rowSize)~=0)
        pairwisePotMat(1,4) = -1.*log(horizontalAMat(currentLabels(patchID+1),currentLabels(patchID)));
    else
        pairwisePotMat(1,4) = 0;    % nothing to the right
    end
    
    
    f2 = sum(pairwisePotMat,2);