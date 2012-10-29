function coefMat = sampleFromMRF(initLabels,Dictionary,rowSize,...
                    colSize,verticalAMat,horizontalAMat,sigma)

% Inputs:
% initLabels - initial labels for each patch. num cols = num patches.
% Dictionary
% rowSize
% colSize
% verticalAMat
% horizontalAMat

%% Init
totPatches = size(initLabels,2);
coefMat = sparse(size(Dictionary,2),totPatches);
listOfNotSampledPatchIndices = 1:totPatches;

%% Sampling
for i = 1:totPatches    
    % pick a random patch which is not yet sampled
    ind = ceil(rand(1)*length(listOfNotSampledPatchIndices));
    currentPatchInd = listOfNotSampledPatchIndices(ind);

    % pick a label for this patch according to the conditional prob distr 
    % based on the neighborhood
    wordIndForCurrPatch = getWordForThisPatch(patchID,Dictionary,verticalAMat,horizontalAMat,...
                    currentLabels,inputData,rowSize,colSize,sigma);
    coefMat(wordIndForCurrPatch,currentPatchInd) = 1;
    
    % mark this patch as sampled. i.e. remove it from the list to be
    % sampled
    % b = a(a~=3);
    listOfNotSampledPatchIndices = ...
        listOfNotSampledPatchIndices(listOfNotSampledPatchIndices~=currentPatchInd);
    remainingPatches = length(listOfNotSampledPatchIndices);
    display('Remaining Patches: ');
    display(remainingPatches);

end
