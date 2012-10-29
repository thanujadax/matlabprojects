function wordInd = getWordForThisPatch(patchID,Dictionary,verticalAMat,horizontalAMat,...
    currentLabels,inputData,rowSize,colSize,sigma)

% inputData - row vector containing the current patch. or 0 if no input
% data is provided

z = getZ(Dictionary,inputData,verticalAMat,horizontalAMat,currentLabels,rowSize,colSize);

% f2 column vector containing pairwise potentials for this patch
% each elements assumes this patch is labeled Dictionary(i)
f2 = getPairwisePotentials(patchID,currentLabels,verticalAMat,horizontalAMat,rowsize);
if(inputData==0)
    % no unary potentials
    condProbs = exp(- f2) ./ z;
else
    % with unary potentials
    f1 = getUnaryPotentialVec(inputData,Dictionary,sigma);
    condProbs = exp(-f1-f2) ./ z;
end

%% sampling from the conditional prob distribution
cumulativeProb = 0;
x = rand(1);
for i = 1:size(Dictionary,2)
    cumulativeProb = cumulativeProb + condProbs(i);
    if(cumulativeProb>=x)
        wordInd = i;
        break;
    end
end
