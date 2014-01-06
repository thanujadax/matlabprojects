function regionRewards = getRegionOverlapScore(ws,labelImg_indexed)

% Inputs:
%   ws
%   labelImg_indexed: each cluster is given a unique integer. off regions
%   are give 0 labe.

% Output:
%   regionRewards: vector containing the incremental active proportion over
%   the inactive proportion as a percentage of the total area.
%   off labels have a score close to -1 and on labels have a score close to
%   +1

numRegions = max(max(ws));

regionRewards = zeros(numRegions,1);

for i=1:numRegions
    regionRewards(i) = getOverlapOn_Off_fraction(i,ws,labelImg_indexed);
end

function regionReward = getOverlapOn_Off_fraction(i,ws,labelImg_indexed)

regionPixInds_logical = (ws==i);
numPixels = sum(sum(regionPixInds_logical));
labelsForRegionPix = labelImg_indexed(regionPixInds_logical);

numActivePix = sum(labelsForRegionPix>0);
numInactivePix = numPixels -numActivePix;

regionReward = (numActivePix - numInactivePix)/numPixels;