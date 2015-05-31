function edgeProbabilities = getEdgeProbabilitiesFromRFC...
            (forestEdgeProb,rawImage,OFR,edgepixels,edgePriors,...
            boundaryEdgeIDs,edgeListInds,numTrees,psuedoEdgeIDs,psuedoEdges2nodes)
        
        

fm = getEdgeFeatureMat(rawImage,edgepixels,OFR,edgePriors,boundaryEdgeIDs,...
    edgeListInds,psuedoEdgeIDs,psuedoEdges2nodes);
[y_h,v] = classRF_predict(double(fm), forestEdgeProb);

edgeProbabilities = v(:,2)./numTrees;
