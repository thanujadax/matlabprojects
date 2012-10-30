function [coefMat,labelVector] = sampleFromMRF(currentLabels,inputData,Dictionary,rowSize,...
                    verticalAMat,horizontalAMat,sigma,lambda)

% Inputs:
% initLabels - initial labels for each patch. num cols = num patches.
% inputData - im2col matrix or vector of zeros if no input is given.
% Dictionary
% rowSize
% colSize
% verticalAMat
% horizontalAMat
% sigma - normalization factor for pixel error. (not used anymore)
% lambda - weighting parameter for unary potential

%% Init
totPatches = size(currentLabels,2);
coefMat = sparse(size(Dictionary,2),totPatches);
listOfNotSampledPatchIndices = 1:totPatches;
labelVector = zeros(1,size(currentLabels,2));

%% Sampling
for i = 1:totPatches    
    % pick a random patch which is not yet sampled
    ind = ceil(rand(1)*length(listOfNotSampledPatchIndices));
    currentPatchInd = listOfNotSampledPatchIndices(ind);

    % pick a label for this patch according to the conditional prob distr 
    % based on the neighborhood
    wordIndForCurrPatch = getWordForThisPatch(currentPatchInd,Dictionary,verticalAMat,horizontalAMat,...
                    currentLabels,inputData(:,currentPatchInd),rowSize,sigma,lambda);
    coefMat(wordIndForCurrPatch,currentPatchInd) = 1;
    labelVector(1,currentPatchInd) = wordIndForCurrPatch; 
    
    % mark this patch as sampled. i.e. remove it from the list to be
    % sampled
    % b = a(a~=3);
    listOfNotSampledPatchIndices = ...
        listOfNotSampledPatchIndices(listOfNotSampledPatchIndices~=currentPatchInd);
    remainingPatches = length(listOfNotSampledPatchIndices);
    display('Remaining Patches: ');
    display(remainingPatches);

end
