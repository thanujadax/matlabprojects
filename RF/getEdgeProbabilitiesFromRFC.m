function edgeProbabilities = getEdgeProbabilitiesFromRFC...
            (forestEdgeProb,rawImage,OFR,edgepixels,edgePriors,...
            boundaryEdgeIDs,edgeListInds)
        
        

fm = getEdgeFeatureMat(rawImage,edgepixels,OFR,edgePriors,boundaryEdgeIDs,edgeListInds);
[y_h,v] = classRF_predict(double(fm), forestEdgeProb);

edgeProbabilities = v(:,2);