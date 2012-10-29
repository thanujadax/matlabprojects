function wordInd = getWordForThisPatch(patchID,Dictionary,verticalAMat,horizontalAMat,...
    currentLabels,inputData,rowSize,colSize,sigma)
% patchID - index number of the current patch according to im2col
% inputData - row vector containing the current patch. or 0 if no input
% data is provided
wordInd = 0;    % initialize
%z = getZ(Dictionary,inputData,verticalAMat,horizontalAMat,currentLabels,rowSize,colSize);

% f2 column vector containing pairwise potentials for this patch
% each elements assumes this patch is labeled Dictionary(i)
dictionarySize = size(Dictionary,2);
f2 = getPairwisePotentials(patchID,currentLabels,dictionarySize,verticalAMat,horizontalAMat,rowSize);
if(inputData==0)
    % no unary potentials
    %condProbs = exp(- f2) ./ z;
    condProbs = exp(- f2);
else
    % with unary potentials
    f1 = getUnaryPotentialVec(inputData,Dictionary,sigma);
    %condProbs = exp(-f1-f2) ./ z;
    condProbs = exp(-f1-f2);
end
z = sum(condProbs);
condProbs = condProbs./z;
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

if(wordInd==0)
    display('ERROR: word not assigned to this block!!')
end
