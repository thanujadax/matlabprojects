function orientationScoreSpace3D = convolveOrientedGauss_P(img,barLength,barWidth,orientations...
            ,sigX,sigY)

% inputs: 
%   img - contains the pixels that should be used for calculating the hough
%   space
%   thresh - only the pixels above this values would be taken into account

% output:
% orientationScoreSpace3D

% each point in concern, votes for all the bars that it is contributing, at once
[numRows numCols] = size(img);
numOrientations = size(orientations,2);

% create the gaussian kernel
bar = zeros(barWidth,barLength);
centerR = (barWidth+1)/2;
centerC = (barLength+1)/2;

gaussBar = gauss2d(bar,[sigX sigY],[centerC centerR]);
% init
orientationScoreSpace3D = zeros(numRows,numCols,numOrientations);

parfor i=1:numOrientations
    img_i = img;
    % get the gaussian kernel
    gaussBar_i = gaussBar;
    % rotate appropriately
    orientation = orientations(i);
    if(orientation>0)
        orientation = 180 - orientation;
        gaussBar_i = imrotate(gaussBar_i,orientation);
    end
    % convolve to get the orientation score map for this orientation i
    convResult = conv2(img_i,gaussBar_i);
    [convRows convCols] = size(convResult);
    rowMargin = (convRows-numRows)/2;
    colMargin = (convCols-numCols)/2;
    startR = rowMargin+1;
    stopR = convRows-rowMargin;
    startC = colMargin+1;
    stopC = convCols-colMargin;
    orientationScoreSpace3D(:,:,i) = convResult(startR:stopR,startC:stopC);
  % progressbar(j/totPoints); % update progress bar
end
orientationScoreSpace3D = orientationScoreSpace3D./...
        (max(max(max(orientationScoreSpace3D))));