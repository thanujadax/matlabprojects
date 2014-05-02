function c_skeleton = getSkeleton(neuron2Dsections)

% Input
%   neuron2Dsections:
%       calculated from the volume segmentation or ground truth

% Output:
%   skeleton
%       each cell corresponds to a different neuron
%       a cell contains the skeletal points of a neuron ordered from section 
%       1 to n in z. In the same section it's ordered according to x and y 
%       coordinates for branched processes


centroid = getCentroid(neuron2Dsections(i))
