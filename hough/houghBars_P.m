function houghSpace3D = houghBars_P(img,barLength,barWidth,orientations,slidingDist)

% parallelized version

% we only consider non-zero (above a certain threshold) pixels

% parameters
grayThresh = 0.2;

[numRows numCols] = size(img);
[pxlR, pxlC] = find(img>grayThresh);        % contains the indicies of the nonzero pixels 
numOrientations = size(orientations,2);
houghSpace3D = zeros(numRows,numCols,numOrientations);

rowMargin = ceil(barLength-1/2)+2;   % allowing a margin in the image
colMargin = ceil(barLength-1/2)+2;

% first implementation for 4 orientations each 45 degrees apart (default)

%pixWeight=zeros(1,numOrientations);


    %orientation = orientations(orientationInd);
for r=rowMargin:slidingDist:numRows-rowMargin
    for c=colMargin:slidingDist:numCols-colMargin
        parfor i=1:numOrientations
        % for the current center (r,c) get the pixel count
        pixWeight = getPixWeightInBar(img,r,c,orientations(i),barLength,barWidth);
        houghSpace3D(r,c,i) = pixWeight;
        end
            
    end
end
   