function twoRegionEdgeLIDs = getRegionsForEdgeLIDs(setOfRegions_eLIDs,numEdges)

% Input:
%   setOfRegions_edgeLIDs - contains edgeLIDs for each region in row vectors

% Output:
%   twoRegionEdgeLIDs - for each edgeLID give the IDs of the two regions it bounds
%       region 0 will be the image border

twoRegionEdgeLIDs = zeros(numEdges,2);

% numRegions = size(setOfRegions_eLIDs,1);


for i=1:numEdges
    % for each edgeLID, look for the regionIDs that it is being used
    [rIDs,~] = find(setOfRegions_eLIDs==i);
    if(~isempty(rIDs))
        numR_forEdge = numel(rIDs);
        % update twoRegionEdgeLIDs
        if(numR_forEdge==2)
            twoRegionEdgeLIDs(i,1:2) = rIDs;
        elseif(numR_forEdge==1)
            twoRegionEdgeLIDs(i,1) = rIDs;
        else
            disp('ERROR1: getRegionsForEdgeLIDs.m')            
        end        
    end
end

% for i=1:numRegions
%     % for each region get the list of edgeLIDs bounding it
%     edgeLIDs_region = setOfRegions_eLIDs(i,:);
%     edgeLIDs_region = edgeLIDs_region(edgeLIDs_region>0);
%     % update the count in twoRegionEdgeLIDs
%     
% end