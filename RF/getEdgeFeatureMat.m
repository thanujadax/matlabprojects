function fm = getEdgeFeatureMat(rawImage,edgepixels,OFR,edgePriors,...
                boundaryEdgeIDs,edgeListInds,psuedoEdgeIDs,psuedoEdges2nodes,...
                membraneProbabilityMap)

% Inputs:
%   rawImage
%   edgepixels
%   OFR - 3D matrix containing the oriented edge filter response

% Output:
%   fm = feature matrix, each row corresponds to an edge. each column
%   corresponds to a feature.
%       

% features:

%   precalculated features:
%   1: precalculated edge prior (input) (1)

%   OFR based features:
%   2-10: |OFR|: mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
%   11-19: +(OFR): mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
%   20-28: -(OFR): mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
%   29-37: OFR: mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
%   38-45: |OFR| sorted normalized histogram (8)

%   raw pixel value based features:
%   46-50: mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (5)
%   51-58: unsorted normalized histogram (8)
%   59: 1 if the edge is a boundary edge. else 0.

% Tot number of features = 59

numFeatures = 68;
numEdges = size(edgepixels,1);
[sizeR,sizeC,numOFRdim] = size(OFR);
% dimVector = 1:numOFRdim;

fm = zeros(numEdges,numFeatures);

parfor i=1:numEdges
    fm_i = zeros(1,numFeatures);
    edgepixels_i = edgepixels(i,:);
    edgepixels_i = edgepixels_i(edgepixels_i>0);
    if(isempty(edgepixels_i))
        edgepixels_i = getPsuedoEdgePixels...
            (i,psuedoEdgeIDs,psuedoEdges2nodes,edgeListInds);
    end
    
    numEdgePixels = numel(edgepixels_i);
    [r,c] = ind2sub([sizeR sizeC],edgepixels_i);
    % 1. precalculated edge prior
    k = 1;
%     fm(i,k) = edgePriors(i);
    fm_i(k) = edgePriors(i);
    % |OFR|max - |OFR|min (1)
    OFR_i = zeros(numEdgePixels,numOFRdim);
    OFRabs_i = zeros(numEdgePixels,numOFRdim);
    OFRpos_i = zeros(numEdgePixels,numOFRdim);
    OFRneg_i = zeros(numEdgePixels,numOFRdim);
    for j=1:numOFRdim
        z = ones(1,numEdgePixels) * j;
        pixInd3D_j = sub2ind([sizeR sizeC numOFRdim],r,c,z);
        OFR_ij = OFR(pixInd3D_j);
        OFR_i(:,j) = OFR_ij;          
        
        OFRabs_i(:,j) = abs(OFR_ij);
        
        OFRpos_ij = OFR_ij(OFR_ij>0);
        if(~isempty(OFRpos_ij))
            numZeroPads = numEdgePixels - numel(OFRpos_ij);
            OFRpos_ij = padarray(OFRpos_ij,[0 numZeroPads],'post');
            OFRpos_i(:,j) = OFRpos_ij;
        end
        
        OFRneg_ij = OFR_ij(OFR_ij<0);
        if(~isempty(OFRneg_ij))
            numZeroPads = numEdgePixels - numel(OFRneg_ij);
            OFRneg_ij = padarray(OFRneg_ij,[0 numZeroPads],'post');
            OFRneg_i(:,j) = OFRneg_ij;
        end              
    end
    % 2-10. |OFR|: mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    k = k + 1;
    K = k:(k+8);
    [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin] ...
        = getMatStats(OFRabs_i);
    fm_i(K) = [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin];
    % 11-19: +(OFR): mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    k=K(end) + 1;
    K = k:(k+8);
    [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin] ...
        = getMatStats(OFRpos_i);
    fm_i(K) = [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin];
    % 20-28: -(OFR): mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    k=K(end) + 1;
    K = k:(k+8);
    [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin] ...
         = getMatStats(OFRneg_i);
    fm_i(K) = [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin];
    % 29-37: OFR: mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    k=K(end) + 1;
    K = k:(k+8);
    [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin] ...
        = getMatStats(OFR_i);
    fm_i(K) = [mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin];
    % 38-45: |OFR| sorted histogram (8)
    k=K(end) + 1;
    K = k:(k+7);
    ofrhist = sort(hist((sum(OFRabs_i,2)),8));
    fm_i(K) = ofrhist./numEdgePixels;
    % 46-54: pixel intensities - mean,max,min,std,median (5)
    k=K(end) + 1;
    K = k:(k+4);
    rawImage_i = rawImage;
    rawImage_i = invertImage(rawImage_i);
    pixIntensityVector = rawImage_i(edgepixels_i);
    [mean,max,min,std,median] ...
     = getVecStats(pixIntensityVector);
    fm_i(K) = [mean,max,min,std,median];
    % 51-58: unsorted normalized histogram of pixel intensities (8)
    intensityHist = sort(hist(pixIntensityVector,8));
    k=K(end) + 1;
    K = k:(k+7);
    fm_i(K) = intensityHist./numEdgePixels;
    % 59: isBoundaryEdge
    k = K(end) + 1;
    K = k;
    fm_i(K) = isBoundaryEdge(i,boundaryEdgeIDs,edgeListInds)
    % 60-68: membrane probabilities (from low-level RFC)
    k=K(end) + 1;
    K = k:(k+4);
    rawImage_i = membraneProbabilityMap;
    rawImage_i = invertImage(rawImage_i);
    pixIntensityVector = rawImage_i(edgepixels_i);
    [mean,max,min,std,median] ...
     = getVecStats(pixIntensityVector);
    fm_i(K) = [mean,max,min,std,median];
    
    fm(i,:) = fm_i;
end


function b_isOnBoundary = isBoundaryEdge(edgeListInd_i,boundaryEdgeIDs,edgeListInds)

edgeID_i = edgeListInds(edgeListInd_i);
if(sum(ismember(boundaryEdgeIDs,edgeID_i))>0)
    b_isOnBoundary = 1;
else
    b_isOnBoundary = 0;
end