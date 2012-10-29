function f2 = getPairwisePotentialVec(Dictionary,currentLabels,verticalAMat,horizontalAMat,rowSize)

% Inputs:
% Dictionary -
% currentLabels - current labels for all the patches (num of cols = num of patches)
% verticalAMat
% horizontalAMat

N = getNeighborhoods(currentLabels,rowSize,colSize);        % TODO
% i th row of N contains the 4 neighbors for patch i as in currentLables
% N(i) = [top, bottom, left, right]

pairwisePotMat = zeros(size(currentLabels,1),4);

for i = 1:size(currentLabels,2)
    condProbGivenTop = verticalAMat(currentLabels(i),currentLables(i-rowSize));
    condProbGivenBottom = verticalAMat(currentLables(i+rowSize),currentLabels(i));
    condProbGivenLeft = horizontalAMat(currentLabels(i),currentLabels(i-1));
    condProbGivenRight = horizontalAMat(currentLabels(i+1),currentLabels(i));
    
end
