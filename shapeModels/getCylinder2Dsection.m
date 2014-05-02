function [s_p,s_s,d] = getCylinder2Dsection(neuronSection2D,before,after)

% Inputs:
%   neuronSection2D
%   before: section i-1
%   after: section i+1

% Outputs:
%   s_p: coordintates of skeletal point of the input 2D neuron section
%       s_p.x, s_p.y, s_p.z
%   s_s: tangential unit vector at the skeletal point
%       s_s.x, s_s.y, s_s.z
%   d: diameter of the largest circle perpendicular to s_s that can be
%   fitted inside the input 2D neuron section

s_p = calCentroid2Dsection();

s_s = calTangentToSkeleton();

d = calCylinderDiameter();