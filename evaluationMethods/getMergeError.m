function [numWrongMerges,numObj] = getMergeError(trueLabelVec,segLabelVec)

% get how many true labels are there

% for each seglabel get how many true labels exist
% anything more than 1 is a wrong merge

trueLabelsUnique = unique(trueLabelVec);

trueLabelsUnique = trueLabelsUnique(trueLabelsUnique>0);

numObj = numel(trueLabelsUnique);

segLabelsUnique = unique(segLabelVec);
segLabelsUnique = segLabelsUnique(segLabelsUnique>0);

numWrongMerges = 0;

for i = 1:numel(segLabelsUnique)
    sLabelPos_logical = (segLabelVec==segLabelsUnique(i));
    tLabels_for_segLabel = trueLabelVec(sLabelPos_logical);
    tLabels_for_segLabel = tLabels_for_segLabel(tLabels_for_segLabel>0);
    numErrors_i = numel(unique(tLabels_for_segLabel)) - 1;
    if(numErrors_i>0)
       numWrongMerges =  numWrongMerges + numErrors_i;
    end
end