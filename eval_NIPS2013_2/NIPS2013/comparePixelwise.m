function [nPixAccTot,mPixMeasuresPerLabel] = comparePixelwise(mGroundTruth,mAlgoOutput)

if(size(mGroundTruth) ~= size(mAlgoOutput))
    disp('Error - GroundTruth and AlgoOutput are not of the same size');
else
    nPixAccTot = 1-sum(sum(mGroundTruth~=mAlgoOutput))/numel(mAlgoOutput);
    
    vLevels = unique(mGroundTruth);
    nClasses = length(vLevels);
    
    mPixMeasuresPerLabel = zeros(nClasses,4);
    for i=1:nClasses
        mTempGroundTruth = mGroundTruth;
        mTempAlgoOutput = mAlgoOutput;

        mTempGroundTruth(mTempGroundTruth ~= vLevels(i)) = 0;
        mTempGroundTruth(mTempGroundTruth == vLevels(i)) = 1;

        mTempAlgoOutput(mTempAlgoOutput ~= vLevels(i)) = 0;
        mTempAlgoOutput(mTempAlgoOutput == vLevels(i)) = 1;

        m1 = (mTempAlgoOutput==1 & mTempGroundTruth==1);
        m2 = (mTempAlgoOutput==0 & mTempGroundTruth==0);

        nTP = sum(sum(m1));
        nTN = sum(sum(m2));
        nFP = sum(sum(mTempAlgoOutput>mTempGroundTruth));
        nFN = sum(sum(mTempAlgoOutput<mTempGroundTruth));

        mPixMeasuresPerLabel(i,1) = nTP/(nTP+nFP); % Precision / Positive predictive value
        mPixMeasuresPerLabel(i,3) = nTP/(nTP+nFN); % Recall / Sensitivity
        mPixMeasuresPerLabel(i,4) = nTN/(nFP+nTN); % Specificity
        mPixMeasuresPerLabel(i,2) = nTN/(nTN+nFN); % Negative predictive value
    end
end

end