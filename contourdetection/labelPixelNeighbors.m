function edgePixLabels = labelPixelNeighbors(edgePixLabels,pixID,edgeLabel,...
                sizeR,sizeC)
% gives the label pixLabel to the current pixel (pixID) and all its
% neighbors
% Inputs:
%   edgePixLabels: N-by-2 matrix where the first column contains N pixelIDs
%   and the second column contains the corresponding edge labels ( 0 if no
%   label is assigned to that pixel)
%   pixID: index of the current pixel in terms of the original image
%   coordinates (from watershed)
%   edgeLabel: label to be assigned to the pixel and its connected pixels
%   sizeR,sizeC: dimensions of the original image (watershed)

listInd = find(edgePixLabels(:,1)==pixID);
if(isempty(listInd))
    % no such pixel found
    disp('WARNINING: pixel not found. labelPixelNeighbors.m %d',pixID);
else
    % pixel found on the list
    currentLabel = edgePixLabels(listInd,2);
    if(currentLabel~=0)
        % pixel label already set
        if(currentLabel~=edgeLabel)
            disp('ERROR: pixel label conflict! labelPixelNeighbors.m %d',pixID);
        end
    else
        % current label is 0. set pixel label
        edgePixLabels(listInd,2) = edgeLabel;
    
        % deal with neighbors
        % get neighbor pixel indices
        neighbors = getNeighbors(pixID,sizeR,sizeC);
        % recursive call for each neighbor 
        numNeighbors = numel(neighbors);
        if(numNeighbors>0)
           for n=1:numNeighbors
               % if the neighbor is an edge pixel, do the recursive labeling
               neighborListInd = find(edgePixLabels(:,1)==neighbors(n));
               if(~isempty(neighborListInd) && edgePixLabels(neighborListInd,2)==0)
                    edgePixLabels = labelPixelNeighbors(edgePixLabels,neighbors(n),edgeLabel,...
                    sizeR,sizeC);
               end

           end
        end
    end
    
end