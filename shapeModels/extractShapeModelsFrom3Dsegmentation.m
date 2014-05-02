function c_shapeModels = extractShapeModelsFrom3Dsegmentation(inputSegmentation)

% Input: 
%   inputSegmentation: ground truth or volume segmentation for which shape
%   models have to be fitted

% Output
%   c_shapeModels:
%       each cell corresponds to a neuron
%       contains a table of [p_i,s_i,d_i] 


% computes generic cylinders 
%   f_i(p_i,s_i,d_i)
%       p_i: position vector (x,y,z) of skeleton at i 
%       s_i: tangential direction of skeleton at point i based on nearby
%       skeletal points
%       d_i: diameter of the largest circle at i, perpendicular to s_i

neuronProperties = getNeuron