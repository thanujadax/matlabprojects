function [newWS, removedWsIDs] = getCorrectedWSregions(ws,removeEdgeIDs,edges2pixels)
% merges ws regions after removing removeEdgeIDs

% Inputs
%   ws - ws transform of the original image
%   removeEdgeIDs - edges to be removed from ws to merge regions
%   edges2pixels - first col has edgeID. other cols have pixelInds

% Outputs
%   newWS - new ws transform after merging ws region due to removal of
%   edges
%   removedWsIDs - list of removed ws region ids due to merging with larger regions
removedWsIDs = [];
expandedWsIDs = [];
newWS = ws;
edgeIDsAll = edges2pixels(:,1);
edgepixelsAll = edges2pixels;
edgepixelsAll(:,1) = [];
[sizeR,sizeC] = size(ws);
% for each removed edge
for i=1:numel(removeEdgeIDs)
% get the two ws regions on either side
    edgeListInd_i = find(edgeIDsAll==removeEdgeIDs(i));
    edgepixels_i = edgepixelsAll(edgeListInd_i,:);
    edgepixels_i = edgepixels_i(edgepixels_i>0);
    regionIDs = getRegionsForEdgePixels(newWS,edgepixels_i,sizeR,sizeC);
% find the bigger region and assign its wsID to the smaller region except
% if its region 1 = border (takes precedence)
    if(ismember(1,regionIDs))
        % assign id 1 to the other region
        regionIDs = sort(regionIDs);
        newWS(newWS==regionIDs(2)) = 1;
        removedWsIDs = [removedWsIDs regionIDs(2)];
        % TODO: assign the edge pixels also 1
        newWS(edgepixels_i) = 1;
    else
        % if one of the regions is already in the expandedWsIDs list, keep
        % expanding that. i.e. assign that id to the other
        if(ismember(regionIDs(1),expandedWsIDs))
            newWS(newWS==regionIDs(2)) = regionIDs(1);
            removedWsIDs = [removedWsIDs regionIDs(2)];
            % TODO: assign the edge pixels also the ws id
            newWS(edgepixels_i) = regionIDs(1); 
        elseif(ismember(regionIDs(2),expandedWsIDs))
            newWS(newWS==regionIDs(1)) = regionIDs(2);
            removedWsIDs = [removedWsIDs regionIDs(1)];
            % TODO: assign the edge pixels also the ws id
            newWS(edgepixels_i) = regionIDs(2);            
        else
            newWS(newWS==regionIDs(2)) = regionIDs(1);
            removedWsIDs = [removedWsIDs regionIDs(2)];
            % TODO: assign the edge pixels also the ws id
            newWS(edgepixels_i) = regionIDs(1);
            expandedWsIDs = [expandedWsIDs regionIDs(1)];
        end
    end

end