function pixIndRect = getEdgeCoverage_rect(edgePixels,edgeLength,...
            sizeR,sizeC,edgeOrientation,barHalfWidth)

% Output:
%   pixIndRect: pix indices for a rectangular area covering the edge in
%   concern

%% get edge center
numEdgePixels = numel(edgePixels);
centerEdgePixPos = floor((numEdgePixels+1)/2);
centerPixInd = edgePixels(centerEdgePixPos);
[rCenter,cCenter] = ind2sub([sizeR sizeC],centerPixInd);

%% define bar dimensions (zero orientation)
barLength = edgeLength;
barWidth = 1 + barHalfWidth*2;
bar = ones(barWidth,barLength);

%%  position+rotate the bar on image grid and get relevant pixels
rotatedBar = imrotate(bar,(-edgeOrientation));
[rotBarCentR,rotBarCentC] = getRotatedBarCenter(rotatedBar);

% calculate shift
shiftR = rCenter - rotBarCentR;
shiftC = cCenter - rotBarCentC;
% shift bar pix inds
barPixInds = find(rotatedBar);
[barR,barC] = ind2sub([sizeR sizeC],barPixInds);
barR = barR + shiftR;
barC = barC + shiftC;

% retain only the positive valued pixels
barR_pos_logical = barR>0;
barC_pos_logical = barC>0;
% when both are positive
pos_ind_logical = barR_pos_logical & barC_pos_logical;
barR_pos = barR(pos_ind_logical);
barC_pos = barC(pos_ind_logical);

pixIndRect = sub2ind([sizeR sizeC],barR_pos,barC_pos);


%% Supplementary Functions
function [rotBarCentR,rotBarCentC] = getRotatedBarCenter(rotatedBar)
[numR,numC] = size(rotatedBar);
rotBarCentR = floor((numR+1)/2);
rotBarCentC = floor((numC+1)/2);