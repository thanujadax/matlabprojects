function FMs_for_directNeighborSet = readFaceFeatures...
                (fm_cellInteriors,thisSectionID,setOfDirectNeighbors_6)
        
% feature vector for face a = {f(A-B),f(A),f(B)}
%   A - thisCell
%   B - neighbors

% Output
%   FMs_for_directNeighborSet: 6 x numFaceFeatures
%       each row corresponds to a face given in the setOfDirectNeighbors_6

fA = fm_cellInteriors(thisSectionID,:);

fAset = repmat(fA,6,1);

fBset = fm_cellInteriors(setOfDirectNeighbors_6,:);

fA_B = fAset - fBset;

FMs_for_directNeighborSet = [fA_B fAset fBset];