function edgeResp = calculateEdgeResp(pixIndOfr,OFR,...
                orientation,eps_orientation,ofr_stepSize)
            
% Output:
%   edgeResp: unweighted signed edge response around the edge in concern
%   obtained from OFR 3D matrix using the pixels in the proximity of the
%   edge (pixIndOfr), for the given 'orientation' which is one of two
%   possible complementary orientations per each edge.

% Inputs:
%   pixIndOfr: pixels in the proximity of the edge, for which
%   OFR: 3D matrix 
%   orientation: orientation of the edge in concern
%   eps_orientaion: orientation +- eps_orientation is considered when
%   calculating the edge response

% which OFR_orientations (OFR(:,:,orientaion_dim)) should be considered

orientation_min = orientation - eps_orientation;
orientation_max = orientation + eps_orientation;

% when orientation min and max are between 0 and 360
if(orientation_min>=0 && orientation_min<=360 && ...
                orientation_max>=0 && orientation_max<=360)
    dimStart = ceil(orientation_min/ofr_stepSize); % min val is 1
    dimStop = ceil(orientation_max/ofr_stepSize); % max val is numDimOFR
    
    dimSeq = dimStart:dimStop;
    
% when orientation_min is less than 0, but max is between 0 and 360
elseif(orientation_min<0 && orientation_max>=0 && orientation_max<=360)
    orientation_min = 360 + orientation_min;
    % orientationSteps from orientation_min upto 360
    dimStart = ceil(orientation_min/ofr_stepSize);
    dimStop = ceil(360/ofr_stepSize);
    dimSeq1 = dimStart:dimStop;
    
    % orientationSteps from 0 to orientation_max
    dimStop = ceil(orientation_max/ofr_stepSize);
    dimSeq2 = 1:dimStop;
    
    dimSeq = [dimSeq1 dimSeq2];

% when orientation_min is between 0 and 360 but max is more than 360
elseif(orientation_min>=0 && orientation_min<=360 && ...
        orientation_max>360)
% when orientation_min is     
else
    disp('ERROR in calculateEdgeResp.m. out of range dimensions!')
end

% Calculate edge response using the sequence of relevant OFR dimensions obtained
edgeResp = edgeRespForDims(OFR,dimSeq,pixIndOfr);

%% Supplementary functions
function edgeResp = edgeRespForDims(OFR,dimSeq,pixIndOfr)
numDimsToLook = numel(dimSeq);
numEdgePix = numel(pixIndOfr);
resp = 0;
for i=1:numDimsToLook
    ofr_dim = OFR(:,:,dimSeq(i));
    resp_dim = ofr_dim(pixIndOfr);
    resp = resp + sum(sum(resp_dim));
end

% normalize by the number of dimensions and pixels
edgeResp = resp/(numEdgePix*numDimsToLook);