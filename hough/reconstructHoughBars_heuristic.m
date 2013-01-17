function output = reconstructHoughBars_heuristic(peaks3D,orientations,barLength,barWidth)
% reconstructs the original image based on the peaks given for each
% orientation
% Parallel

% Inputs:
%   peaks3D - a 3D array [row col orientation] containing the votes
%   orientations - e.g. [0 45 90 135]. can be any value from 1 to 180.
%   barLength - 
%   barWidth -

% First pick the bars where we have a high confidence. 
% Iterate:
%   Grow the configuration so that it minimizes some energy function
%   Prune the existing configuration to minimize the energy

%% Parameters
thresh = 0.8;           % threshold to pick the most confident points

%% Reconstruction algorithm - init
[numRows numCols numOrientations] = size(peaks3D);
output = zeros(numRows,numCols,3);

% Pick the bars with high confidence
peaks3D = peaks3D./max(max(max(peaks3D))); % normalize
goodPoints = find(peaks3D>thresh);
output = placeBars(output,goodPoints,peaks3D,orientations,barLength,barWidth);



%% plots
figure(700);imagesc(output);colormap(hsv);
figure(701);imagesc(output(:,:,3));title('votes');
figure(702);imagesc(output(:,:,1));title('orientations')

