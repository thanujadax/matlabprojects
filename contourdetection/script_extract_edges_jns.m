% extracting edges and junctions from WS
imIn = imread('stem_256x_t02_V.png');
[sizeR,sizeC] = size(imIn);
ws = watershed(imIn);
%figure;imagesc(ws)

% the edges (watershed boundaries) are labeled 0
% extract those
ind0 = find(ws==0);
% [r0,c0] = ind2sub([sizeR sizeC],ind0);
wsBoundaries = zeros(sizeR,sizeC);
wsBoundaries(ind0) = 1;
figure;imshow(wsBoundaries);
title('watershed boundaries')

%% extracting junctions
% look at the 4 neighborhood of each pixel
fourNH = zeros(size(imIn));
numBoundaryPixels = numel(ind0);
for i=1:numBoundaryPixels
    % calculate n.o. 4 neighbors
    ind = ind0(i);
    [r c] = ind2sub([sizeR sizeC],ind);
    nh = zeros(1,4);
    
    if((r-1)>0)
        nh(1) = wsBoundaries(r-1,c);
    end
    if((r+1)<=sizeR)
        nh(2) = wsBoundaries(r+1,c);
    end
    if((c-1)>0)
        nh(3) = wsBoundaries(r,c-1);
    end
    if((c+1)<=sizeC)
        nh(4) = wsBoundaries(r,c+1);
    end
    
    fourNH(ind) = sum(nh);
end
% get the pixels which are having a 4NH > 2
ind4J = find(fourNH>2);
% visualize junctions
wsJ = zeros(sizeR,sizeC);
wsJ(ind4J) = 1;
wsVis = zeros(sizeR,sizeC,3);
wsVis(:,:,3) = wsBoundaries;
wsVis(:,:,1) = wsJ;
figure;imshow(wsVis);
title('Junctions from WS')



