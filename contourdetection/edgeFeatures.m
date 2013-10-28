function fm = edgeFeatures(rawImage,edgepixels,OFR,edgePrior)

% Inputs:
%   edgepixels - each row corresponds to an edge. gives the pixelInds as a
%   row vector
%   OFR - the 3D matrix containing OFR. 
%   edgePrior - pre-calculated edgePriors if available.

numEdges = size(edgepixels,1);
[sizeR,sizeC,numOFRdim] = size(OFR);

% initialize feature matrix
fm = zeros(numEdges,numOFRdim);

for i=1:numEdges
    % get all edge pixels
    % get OFR for each pixel (all dimensions)
    clear edgePixInds_i
    edgePixInds_i = edgepixels(i,:);
    edgePixInds_i = edgePixInds_i(edgePixInds_i>0);
    clear edgePixInd_r edgePixInd_c 
    [edgePixInd_r,edgePixInd_c] = ind2sub([sizeR sizeC],edgePixInds_i);
    numEdgePix_i = numel(edgePixInds_i);
    
    edgeOFR_i = zeros(numEdgePix_i,numOFRdim);
    for dim=1:numOFRdim
        ofr_dim = zeros(sizeR,sizeC);
        ofr_dim = OFR(:,:,dim);
        edgeOFR_i(:,dim) = ofr_dim(edgePixInds_i);
        fm(i,dim) = max(edgeOFR_i(:,dim));
        k = 1;
        fm(i,(dim+numOFRdim)) = min(edgeOFR_i(:,dim));
        k = k + 1;
        fm(i,(dim+numOFRdim*k)) = median(edgeOFR_i(:,dim));
        k = k + 1;
        fm(i,(dim+numOFRdim*k)) = mean(edgeOFR_i(:,dim));
        k = k + 1;
        fm(i,(dim+numOFRdim*k)) = std(edgeOFR_i(:,dim));
        k = k + 1;
        fm(i,(dim+numOFRdim*k)) = mode(edgeOFR_i(:,dim));
    end
    
    numDimFeatures = numOFRdim * k;
    
    max_ofr_i = max(fm(i,1:numDimFeatures));
    mode_ofr_i = mode(fm(i,1:numDimFeatures));
    min_ofr_i = min(fm(i,1:numDimFeatures));
    diff_ofr_i = max_ofr_i - min_ofr_i;
    
    k = 1;
    fm(i,numDimFeatures+k) = max_ofr_i;
    k = k + 1;
    fm(i,numDimFeatures+k) = min_ofr_i;
    k = k + 1;
    fm(i,numDimFeatures+k) = mode_ofr_i;
    k = k + 1;
    fm(i,numDimFeatures+k) = diff_ofr_i;
    
    k = k + 1;
    fm(i,numDimFeatures+k) = numEdgePix_i;
    
    k = k + 1;
    fm(i,numDimFeatures+k) = edgePrior(i);  % pre-calculated edge prior (input)
    
end

