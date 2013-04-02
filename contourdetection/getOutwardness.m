function outwardness = getOutwardness(theta,alpha,C)
% returns the degree of outwardness for an edge, scaled by parameter C

% Inputs:
%   theta - angle of the edge in concern based on OFR (array)
%   alpha - angle of the edge based on its position on the graph (array)
%   C - scaling factor

dta = alpha - theta;
outwardness = cosd(dta).*C;