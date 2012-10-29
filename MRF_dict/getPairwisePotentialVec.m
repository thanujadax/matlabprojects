function f2 = getPairwisePotentialVec(currentLabels,verticalAMat,horizontalAMat,rowSize)

% Inputs:
% currentLabels - current labels for all the patches (num of cols = num of patches)
% verticalAMat
% horizontalAMat

% Output:
% f2 - pairwise potentials as a col vector for each patch

% N = getNeighborhoods(currentLabels,rowSize,colSize);
% i th row of N contains the 4 neighbors for patch i as in currentLables
% N(i) = [top, bottom, left, right]
% However, in this code, the neighborhood elements are directly accessed
% from currentLabels matrix.

pairwisePotMat = zeros(size(currentLabels,1),4);

for i = 1:size(currentLabels,2)
    % condCostGivenTop
    if(i - rowSize > 0)
        pairwisePotMat(i,1) = -1.*log(verticalAMat(currentLabels(i),currentLables(i-rowSize)));
    else
        pairwisePotMat(i,1) = 0;    % nothing on top
    end
    % condCostGivenBottom
    if(i+rowSize <= size(currentLabels,2))
        pairwisePotMat(i,2) = -1.*log(verticalAMat(currentLables(i+rowSize),currentLabels(i)));
    else
        pairwisePotMat(i,2) = 0;    % nothing to the bottom
    end
    % condCostGivenLeft
    if(i-1 > 0 && mod(i,rowSize)>1)
        pairwisePotMat(i,3) = -1.*log(horizontalAMat(currentLabels(i),currentLabels(i-1)));
    else
        pairwisePotMat(i,3) = 0;    % nothing to the left
    end
    % condCostGivenRight
    if(i+1<=size(currentLabels,2) && mod(i,rowSize)~=0)
        pairwisePotMat(i,4) = -1.*log(horizontalAMat(currentLabels(i+1),currentLabels(i)));
    else
        pairwisePotMat(i,4) = 0;    % nothing to the right
    end
    
end



f2 = sum(pairwisePotMat,2);