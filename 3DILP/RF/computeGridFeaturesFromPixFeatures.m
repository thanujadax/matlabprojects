function fm_gridForSlice_i = computeGridFeaturesFromPixFeatures...
                            (fm_pixCell)
                        
% Compute the feature matrix for the gridCell given by the set of features
% of the relevant pixels

% Inputs:
%   fm_pixCell - 3D matrix containing the features for the grid pixels

[numPixR,numPixC,numPixFeatures] = size(fm_pixCell);
numPix = numPixR * numPixC;

% min,max,mean,var,median of all pixel features
numGridCellFeatures = numPixFeatures * 5;

fm_gridForSlice_i = zeros(1,numGridCellFeatures);
% stopInd = 0;
f_ind = 0;
for i=1:numPixFeatures
%    startInd = stopInd + 1;
%    stopInd = stopInd + 5;

   pixFeature_i = fm_pixCell(:,:,i);
   nanInds = isnan(pixFeature_i);
   if(sum(sum(nanInds))>0)
       pixFeature_i(nanInds) = 0;
   end
   pixFeatureVector_i = reshape(pixFeature_i,[1 numPix]);
   % min
   f_ind = f_ind + 1;
   fm_gridForSlice_i(f_ind) = min(pixFeatureVector_i);
   % max
   f_ind = f_ind + 1;
   fm_gridForSlice_i(f_ind) = max(pixFeatureVector_i);
   % mean
   f_ind = f_ind + 1;
   fm_gridForSlice_i(f_ind) = mean(pixFeatureVector_i);
   % var
   f_ind = f_ind + 1;
   fm_gridForSlice_i(f_ind) = var(pixFeatureVector_i);
   % median
   f_ind = f_ind + 1;
   fm_gridForSlice_i(f_ind) = median(pixFeatureVector_i);
   
    
end