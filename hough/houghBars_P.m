function houghSpace3D = houghBars_P(img,barLength,barWidth,orientations,slidingDist)

% parallelized version

% we only consider non-zero (above a certain threshold) pixels

% parameters
%grayThresh = 0.2;

[numRows numCols] = size(img);
%[pxlR, pxlC] = find(img>grayThresh);        % contains the indicies of the nonzero pixels 
numOrientations = size(orientations,2);
houghSpace3D = zeros(numRows,numCols,numOrientations);

rowMargin = ceil(barLength-1/2)+2;   % allowing a margin in the image
colMargin = ceil(barLength-1/2)+2;

% first implementation for 4 orientations each 45 degrees apart (default)

%pixWeight=zeros(1,numOrientations);

totIter = (numRows - 2*rowMargin)*numOrientations; % for the progress bar
counter = 0; % for the progress bar
progressbar('Calculating 3D Hough space'); % Create figure and set starting time

    %orientation = orientations(orientationInd);
display('houghBars_P: parallel computation of the Hough space')
for i=1:numOrientations
    for r=rowMargin:numRows-rowMargin
        orientation = orientations(i);
        parfor c=colMargin:numCols-colMargin
            % for the current center (r,c) get the pixel count
            pixWeight = getPixWeightInBar(img,r,c,orientation,barLength,barWidth);
            houghSpace3D(r,c,i) = pixWeight;
        end
        counter = counter + 1;    
        progressbar(counter/totIter); % Update progress bar
    end
end

%progressbar(1);
   