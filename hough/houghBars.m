function houghSpace3D = houghBars(img,barLength,barWidth,orientations,slidingDist)

% TODO: parallelize

% we only consider non-zero (above a certain threshold) pixels

% parameters
grayThresh = 0.3;

[numRows numCols] = size(img);
[pxlR, pxlC] = find(img>grayThresh);        % contains the indicies of the nonzero pixels 
numOrientations = size(orientations,2);
houghSpace3D = zeros(numRows,numCols,numOrientations);

rowMargin = ceil(barLength-1/2)+2;   % allowing a margin in the image
colMargin = ceil(barLength-1/2)+2;

% first implementation for 4 orientations each 45 degrees apart (default)

totIter = (numRows - 2*rowMargin)*numOrientations; % for the progress bar
counter = 0; % for the progress bar
progressbar('Calculating 3D Hough space') % Create figure and set starting time

for orientationInd=1:numOrientations
    orientation = orientations(orientationInd);
    for r=rowMargin:slidingDist:numRows-rowMargin
        for c=colMargin:slidingDist:numCols-colMargin
            % for the current center (r,c) get the pixel count
            pixWeight = getPixWeightInBar(img,r,c,orientation,barLength,barWidth);
            houghSpace3D(r,c,orientationInd) = pixWeight;
            
        end
    end
    counter = counter + 1;    
    progressbar(counter/totIter) % Update progress bar 
end
progressbar(1);