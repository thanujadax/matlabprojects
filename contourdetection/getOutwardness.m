function outwardness = getOutwardness(theta,alpha)
% returns the degree of outwardness for an edge

% Inputs:
%   theta - angle of the edge in concern based on OFR (array)
%   alpha - angle of the edge based on its position on the graph (array)


dta = alpha - theta;
outwardness = cosd(dta);