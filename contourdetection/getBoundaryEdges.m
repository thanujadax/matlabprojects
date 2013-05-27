function boundaryEdges = getBoundaryEdges(ofr,barLength,ws_edgePixels)

% Inputs:
%   ofr - absolute value of OFR
%   barLength - 
%   wsEdgePixels -

visualize = 1;

[sizeR,sizeC] = size(ofr);
marginInd = floor(barLength/2);

% get all possible boundary pixels from OFR
% separately for each boundary
topPixInd_col = find(ofr(marginInd,:)>0);
numTopPix = numel(topPixInd_col);
topPixInd_row = ones(numTopPix,1) .* marginInd;
topPixInd = sub2ind([sizeR sizeC],topPixInd_row,topPixInd_col');

botPixInd_col = find(ofr((sizeR-marginInd),:)>0);
numBotPix = numel(botPixInd_col);
botPixInd_row = ones(numBotPix,1) .* (sizeR-marginInd);
botPixInd = sub2ind([sizeR sizeC],botPixInd_row,botPixInd_col');

leftPixInd_row = find(ofr(:,marginInd)>0);
numLeftPix = numel(leftPixInd_row);
leftPixInd_col = ones(numLeftPix,1) .* marginInd;
leftPixInd = sub2ind([sizeR sizeC],leftPixInd_row,leftPixInd_col);

rightPixInd_row = find(ofr(:,(sizeC-marginInd)));
numRightPix = numel(rightPixInd_row);
rightPixInd_col = ones(numRightPix,1) .* (sizeC-marginInd);
rightPixInd = sub2ind([sizeR sizeC],rightPixInd_row,rightPixInd_col);

% out of the all possible boundary pixels from OFR, extract the ones that
% correspond to WS edges
edgeListInd_top = getEdgeListIndForPixInds(topPixInd,ws_edgePixels);
edgeListInd_bot = getEdgeListIndForPixInds(botPixInd,ws_edgePixels);
edgeListInd_left = getEdgeListIndForPixInds(leftPixInd,ws_edgePixels);
edgeListInd_right = getEdgeListIndForPixInds(rightPixInd,ws_edgePixels);

if(visualize)
    % visualize detected boundary edges
    imgTmp = zeros(sizeR,sizeC);
    % color all edge pixels : 1
    edgepix_all = ws_edgePixels(ws_edgePixels>0);
    imgTmp(edgepix_all) = 1;
    % color top and bottom edges
    % top 0.3
    numTopEdges = numel(edgeListInd_top);
    for i=1:numTopEdges
        clear edgepix_i
        edgepix_i = ws_edgePixels(edgeListInd_top(i),:);
        edgepix_i = edgepix_i(edgepix_i>0);
        imgTmp(edgepix_i) = 0.3;
    end
    % bottom 0.3
    numBotEdges = numel(edgeListInd_bot);
    for i=1:numBotEdges
        clear edgepix_i
        edgepix_i = ws_edgePixels(edgeListInd_bot(i),:);
        edgepix_i = edgepix_i(edgepix_i>0);
        imgTmp(edgepix_i) = 0.3;
    end
    % color left and right edges
    % left 0.6
    numLeftEdges = numel(edgeListInd_left);
    for i=1:numLeftEdges
        clear edgepix_i
        edgepix_i = ws_edgePixels(edgeListInd_left(i),:);
        edgepix_i = edgepix_i(edgepix_i>0);
        imgTmp(edgepix_i) = 0.6;
    end
    % right 0.6
    numRightEdges = numel(edgeListInd_right);
    for i=1:numRightEdges
        clear edgepix_i
        edgepix_i = ws_edgePixels(edgeListInd_right(i),:);
        edgepix_i = edgepix_i(edgepix_i>0);
        imgTmp(edgepix_i) = 0.6;
    end
    figure;imagesc(imgTmp);title('boundary edges');
end
boundaryEdges = 0;