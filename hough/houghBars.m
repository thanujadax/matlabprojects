function houghSpace3D = houghBars(img,barLength,barWidth,orientations,slidingDist)

% TODO: parallelize

% we only consider non-zero (above a certain threshold) pixels

% parameters
grayThresh = 0.3;

[numRows numCols] = size(img);
[pxlR, pxlC] = find(img>grayThresh);        % contains the indicies of the nonzero pixels 
numOrientations = size(orientations,2);
houghSpace3D = zeros(numRows,numCols,numOrientations);

rowMargin = ceil(barLength/2);   % allowing a margin in the image
colMargin = ceil(barLength/2);

% first implementation for 4 orientations each 45 degrees apart (default)

for orientationInd=1:numOrientations
    for r=rowMargin:slidingDist:numRows-rowMargin
        for c=colMargin:slidingDist:numCols-colMargin
            % for the current center (r,c) get the pixel count
            pixelCount = pixelCountInBar(img,r,c,orientations(orientationInd)...
                    ,barLength,barWidth,grayThresh);
            houghSpace3D(r,c,orientationInd) = pixelCount;
            
        end
    end
    
end