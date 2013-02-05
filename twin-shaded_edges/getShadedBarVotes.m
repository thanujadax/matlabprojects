function houghSpace3D = getShadedBarVotes(img,barLength,barWidth,orientations,offWidth,slidingDist)

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

totIter = (numRows - 2*rowMargin)*numOrientations; % for the progress bar
counter = 0; % for the progress bar
progressbar('Calculating 3D vote space'); % Create figure and set starting time

    %orientation = orientations(orientationInd);
display('getShadedBarVotes: parallel computation of the 3D vote space')
for i=1:numOrientations
    for r=rowMargin:numRows-rowMargin
        orientation = orientations(i);
        parfor c=colMargin:numCols-colMargin
            % for the current center (r,c) get the pixel count
            pixWeight = getTwinShadedBarWeight(img,r,c,orientation,barLength,barWidth);
            %pixWeight = get010BarWeight(img,r,c,orientation,barLength,barWidth,offWidth);
            houghSpace3D(r,c,i) = pixWeight;
        end
        counter = counter + 1;    
        progressbar(counter/totIter); % Update progress bar
    end
end



progressbar(1);