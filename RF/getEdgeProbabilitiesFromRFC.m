function edgeProbabilities = getEdgeProbabilitiesFromRFC...
            (forestEdgeProb,rawImage,OFR,edgepixels,edgePriors)
        
        

fm = getEdgeFeatureMat(rawImage,edgepixels,OFR,edgePriors);
[y_h,v] = classRF_predict(double(fm), forestEdgeProb);

edgeProbabilities = v(:,2);