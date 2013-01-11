function voteMat = peakProcessing(voteMat,orientation,barLength,barWidth,method)
% works independently for each orientation
% for the given orientation, decide which peaks to plot

% Inputs:
%   voteMat - 2D vote matrix for the given orientation 
%   orientation - in degrees. 0 to 180
%   barLength
%   barWidth
%   method - NMS 