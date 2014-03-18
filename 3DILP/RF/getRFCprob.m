function predictedProbs = getRFCprob...
                                (fm,forestMat,numTrees)

% Inputs:
%   fm - 3D feature matrix. Each entry in the 3rd dim corresponds to a
%   feature
%   forest - trained RFC
%   numTrees -

% Output:
%   gridCellInteriorUnary - the probabily that each gridCell is neuron
%   interior

% read precomputed featureMats for each image


% get features for gridCells

% Run classifier
[y_h,v] = classRF_predict(double(fm), forestMat);

predictedProbs = v(:,2)./numTrees;