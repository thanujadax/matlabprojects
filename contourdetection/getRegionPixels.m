function regionPixels = getRegionPixels(regionList,wsIDsForRegions,ws)

% Inputs:
%   regionList: list of regionIDs
%   wsIDsForRegions
%   ws

% Output:
%   regionPixels: internal pixels of all the regions in the list

regionPixels = [];

numRegions = numel(regionList);

for i=1:numRegions
    % get internal pixel of region i
    wsID = wsIDsForRegions(regionList(i));
    internalPix_i = find(ws==wsID);
    % append to list of internal pixels
    regionPixels = [regionPixels; internalPix_i];
end