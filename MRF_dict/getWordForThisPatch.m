function wordInd = getWordForThisPatch(patchID,Dictionary,verticalAMat,horizontalAMat,...
    currentLabels,inputData,rowSize,colSize)

z = getZ(Dictionary,inputData,verticalAMat,horizontalAMat,currentLabels,rowSize,colSize);

f2 = getPairwisePotentials(patchID,currentLabels,verticalAMat,horizontalAMat,rowsize);
if(inputData==0)
    % no unary potentials
    condProbs = exp(- f2) ./ z;
else
    % with unary potentials
    condProbs = exp(-f1-f2) ./ z;
end