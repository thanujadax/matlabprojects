% extracting edges and junctions from WS
imIn = imread('stem_256x_t02_V.png');
ws = watershed(imIn);
%figure;imagesc(ws)

% the edges (watershed boundaries) are labeled 0
% extract those
ind0 = find(ws==0);
wsBoundaries = zeros(size(imIn));
wsBoundaries(ind0) = 1;
figure;imshow(wsBoundaries);

% calculate how many neighbors each boundary pixel has
% get the indices of the boundary pixels
% sum the pixels of a 3x3 neighborhood
convFilter = ones(3);
convOut = conv2(wsBoundaries,convFilter);
figure;imagesc(convOut);

