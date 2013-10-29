function fm = getEdgeFeatureMat(rawImage,edgepixels,OFR,edgePriors)

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

% Tot number of features = 58

numFeatures = 58;
numEdges = size(edgepixels,1);
[sizeR,sizeC,numOFRdim] = size(OFR);
dimVector = 1:numOFRdim;

fm = zeros(numEdges,numFeatures);

for i=1:numEdges
    edgepixels_i = edgepixels(i,:);
    edgepixels_i = edgepixels_i(edgepixels_i>0);
    numEdgePixels = numel(edgepixels_i);
    % 1. precalculated edge prior
    k = 1;
    fm(i,k) = edgePriors(i);      
    % |OFR|max - |OFR|min (1)
    OFR_i = zeros(numEdgePixels,numOFRdim);
    OFRabs_i = zeros(numEdgePixels,numOFRdim);
    OFRpos_i = zeros(numEdgePixels,numOFRdim);
    OFRneg_i = zeros(numEdgePixels,numOFRdim);
    for j=1:numOFRdim
        z = ones * j;
        pixInd3D_j = sub2ind([sizeR sizeC numOFRdim],r,c,z);
        OFR_ij = OFR(pixInd3D_j);
        OFR_i(:,j) = OFR_ij;          
        
        OFRabs_i(:,j) = abs(OFR_ij);
        
        OFRpos_ij = OFR_ij(OFR_ij>0);
        numZeroPads = numEdgePixels - numel(OFRpos_ij);
        OFRpos_ij = padarray(OFRpos_ij,[numZeroPads 0],'post');
        OFRpos_i(:,j) = OFRpos_ij;
        
        OFRneg_ij = OFR_ij(OFR_ij<0);
        numZeroPads = numEdgePixels - numel(OFRneg_ij);
        OFRneg_ij = padarray(OFRneg_ij,[numZeroPads 0],'post');
        OFRneg_i(:,j) = OFRneg_ij;
              
    end
    % 2-10. |OFR|: mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    K = k:(k+8);
    fm(i,K) = getMatStats(OFRabs_i);
    % 11-19: +(OFR): mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    k=K(end) + 1;
    K = k:(k+8);
    fm(i,K) = getMatStats(OFRpos_i);
    % 20-28: -(OFR): mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    k=K(end) + 1;
    K = k:(k+8);
    fm(i,K) = getMatStats(OFRneg_i);
    % 29-37: OFR: mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (9)
    k=K(end) + 1;
    K = k:(k+8);
    fm(i,K) = getMatStats(OFR_i);
    % 38-45: |OFR| sorted histogram (8)
    k=K(end) + 1;
    K = k:(k+7);
    ofrhist = sort(hist(OFRabs_i,8));
    fm(i,K) = ofrhist./sum(ofrhist);
    % 46-54: pixel intensities - mean,max,min,stdMean,stdMax,stdMin,medianMean,medianMax,medianMin (5)
    k=K(end) + 1;
    K = k:(k+4);
    pixIntensityVector = rawImage(edgepixels_i);
    fm(i,K) = getVecStats(pixIntensityVector);
    % 51-58: unsorted normalized histogram of pixel intensities (8)
    intensityHist = sort(hist(pixIntensityVector,8));
    k=K(end) + 1;
    K = k:(k+7);
    fm(i,K) = intensityHist./sum(intensityHist);
end
