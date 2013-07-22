function offEdgeListIDs = getUnOrientedEdgeIDs(edgepixels,...
                lenThresh,orientationScoresMax,sizeR,sizeC)
% returns the IDs of edges that have sharp changes of orientation (~180)
% along its pixels

% Inputs:
%   

% get the list of edges with the corresponding set of pixels for each edge
% for each edge of each pixel, get the maxRespOrientation
% get the diff of the maxRespOrientation values
% if it has a diff around 180 degrees, discard that edge

%% parameters
upThresh = 260;  % orientation difference upper bound
downThresh = 100;% orientation difference lower bound
alphaDiffMaxThresh = 1;  % when the difference of alphas for the 2 edges are
% very small, a stricter range of thresholds are used to improve accuracy
upThreshLowAlpha = 285;
downThreshLowAlpha = 75;
% obvious misorientations (close to 180) - including pixels close to nodes
upThreshLowAlpha_super = 225;
downThreshLowAlpha_super = 135;
%%
offEdgeListIDs = [];
[numEdges,maxPixelsPerEdge] = size(edgepixels);
edgeLengths = zeros(numEdges,1);
for i = 1:numEdges
    edgeLengths(i) = sum((edgepixels(i,:))>0);
end
edgeListIndToExamine = find(edgeLengths<lenThresh);
k = 0;
for i = 1:numel(edgeListIndToExamine)
    clear edgePixels_i;
    edgeListID = edgeListIndToExamine(i); % row number for the edge in edges2pixels
    if(edgeListID==541)
        abc = 1;
    end
    edgePixels_i = edgepixels(edgeListID,:); % set of pixels for 
                % this edge, padded with zero entries
    edgePixels_i = edgePixels_i(edgePixels_i>0); % remove the zero padding
    pixelOrientations = orientationScoresMax(edgePixels_i).*360; % orientation values for the pixels
    % check if the edge contains a unacceptable set of orientations
    % if yes, add it to the list of offEdgeIDs    
    pxOri_diff = abs(diff(pixelOrientations));
    pxEnd_diff = abs(pixelOrientations(1)-pixelOrientations(end));
    % calculate the orientation differences on the actual graph edge
    [pixr,pixc] = ind2sub([sizeR,sizeC],edgePixels_i);
    orientations_alpha = atan2d(pixr,pixc);
    pxOri_alpha_diff = abs(diff(orientations_alpha));
    maxAlphaDiff = max(pxOri_alpha_diff);
    numDiffs = numel(pxOri_diff);
    clear misOrientations;
    if(maxAlphaDiff<alphaDiffMaxThresh)
        % first check if there's a deviation close to 180
%         misOrientations = intersect(pxOri_diff(pxOri_diff>downThreshLowAlpha_super),...
%                 pxOri_diff(pxOri_diff<upThreshLowAlpha_super));
        misOrientations = intersect(pxEnd_diff(pxEnd_diff>downThreshLowAlpha_super),...
                pxEnd_diff(pxEnd_diff<upThreshLowAlpha_super));
        if(isempty(misOrientations))   
            % leave out the pixels close to the nodes.
            if(numDiffs>4)
                pxOri_diff(end) = [];
                pxOri_diff(1) = [];
                pxOri_diff(end) = [];
                pxOri_diff(1) = [];
                misOrientations = intersect(pxOri_diff(pxOri_diff>downThreshLowAlpha),...
                    pxOri_diff(pxOri_diff<upThreshLowAlpha));
            elseif(numDiffs>2)
                pxOri_diff(end) = [];
                pxOri_diff(1) = [];
                misOrientations = intersect(pxOri_diff(pxOri_diff>downThreshLowAlpha),...
                    pxOri_diff(pxOri_diff<upThreshLowAlpha));
            else
                % edge has less than 4 pixels. can't apply the above exclusion condition.
                % misOrientations = [];
                misOrientations = intersect(pxOri_diff(pxOri_diff>downThreshLowAlpha),...
                    pxOri_diff(pxOri_diff<upThreshLowAlpha));
            end
        end
    else
        misOrientations = intersect(pxOri_diff(pxOri_diff>downThresh),pxOri_diff(pxOri_diff<upThresh));
    end
    % misOrientations should contain an element if this edge is
    % misoriented.
    misOriented = numel(misOrientations);
    if(misOriented>0)
        % edge_i is misoriented
        % add to the off list
        k = k + 1;
        offEdgeListIDs(k) = edgeListID;
    end
end
