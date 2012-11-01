function f2 = getPairwisePotentials(patchID,currentLabels,dictionarySize,...
                    verticalAMat,horizontalAMat,rowSize)
                
% normalize association matrices if they are not already normalized
verticalAMat_colNormalized = normalizeAMat(verticalAMat,1);
verticalAMat_rowNormalized = normalizeAMat(verticalAMat,0);

horizontalAMat_colNormalized = normalizedAMat(horizontalAMat,1);
horizontalAMat_rowNormalized = normalizedAMat(horizontalAMat,0);

% for a single patch considering its neighborhoot
pairwisePotMat = zeros(dictionarySize,4);
for i = 1:dictionarySize

    if(patchID - rowSize > 0)
        % pairwisePotMat(i,1) = -1.*log(verticalAMat(i,currentLabels(patchID-rowSize)));
        pairwisePotMat(i,1) = -1.*log(verticalAMat_colNormalized(i,currentLabels(patchID-rowSize)));
    else
        pairwisePotMat(i,1) = 0;    % nothing on top
    end
    % condCostGivenBottom
    if(patchID+rowSize <= size(currentLabels,2))
        pairwisePotMat(i,2) = -1.*log(verticalAMat_rowNormalized(currentLabels(patchID+rowSize),i));
    else
        pairwisePotMat(i,2) = 0;    % nothing to the bottom
    end
    % condCostGivenLeft
    if(patchID-1 > 0 && mod(patchID,rowSize)>1)
        pairwisePotMat(i,3) = -1.*log(horizontalAMat_colNormalized(i,currentLabels(patchID-1)));
    else
        pairwisePotMat(i,3) = 0;    % nothing to the left
    end
    % condCostGivenRight
    if(patchID+1<=size(currentLabels,2) && mod(patchID,rowSize)~=0)
        pairwisePotMat(i,4) = -1.*log(horizontalAMat_rowNormalized(currentLabels(patchID+1),i));
    else
        pairwisePotMat(i,4) = 0;    % nothing to the right
    end
    
end
f2 = sum(pairwisePotMat,2);