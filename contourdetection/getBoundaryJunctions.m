function boundaryJunctions = getBoundaryJunctions(edgepixels,junctionPixels,margin,...
    sizeR,sizeC)
boundaryJunctions = [];
tmp = zeros(sizeR,sizeC);
% get the edges placed in the margins
topR = 1:margin;
allCol = 1:sizeC;
bottomR = (sizeR-margin+1):sizeR;
leftC = 1:margin;
rightC = (sizeC-margin+1):sizeC;
%% get edges - top margin
tmp(topR,:) = 1;
topMarginPixInds = find(tmp);
edgelist_top = [];
for i= 1:numel(topMarginPixInds)
    clear edgelist_tmp
    [edgelist_tmp,~] = find(edgepixels==topMarginPixInds(i));
    edgelist_top = [edgelist_top; edgelist_tmp];
end
edgelist_top = unique(edgelist_top);
% which of these edges are parallel to the boundary
parallelEdgeList_top = [];
for i = 1:numel(edgelist_top)
    pixind = edgepixels(edgelist_top(i),:);
    clear pixr
    clear pixc
    [pixr,pixc] = ind2sub([sizeR sizeC],pixind); 
    pixr = pixr(pixr>0);
    % how many pixels are on the same row
    clear m
    clear s
    m = mode(pixr);
    s = sum(pixr==m);
    if(s>2)
        parallelEdgeList_top = [parallelEdgeList_top;edgelist_top(i)];
    end
end
% test
tmp = zeros(sizeR,sizeC);
tmp(junctionPixels) = 0.6;
pos_edgepix = edgepixels(edgepixels>0);
tmp(pos_edgepix) = 0.4;
topmarginpix = edgepixels(parallelEdgeList_top,:);
topmarginpix_pos = topmarginpix(topmarginpix>0);
tmp(topmarginpix_pos) = 0.8;
figure;imagesc(tmp);
%% bottom margin
tmp = zeros(sizeR,sizeC);
tmp(bottomR,:) = 1;
botMarginPixInds = find(tmp);
%%  left margin
tmp = zeros(sizeR,sizeC);
tmp(leftC,:) = 1;
leftMarginPixInds = find(tmp);
%%  right margin
tmp = zeros(sizeR,sizeC);
tmp(rightC,:) = 1;
rightMarginPixInds = find(tmp);

%% get which edges are there
% out of those edges, get which edges have parallel components to the
% boundaries
% which junctions do these edges connect to

