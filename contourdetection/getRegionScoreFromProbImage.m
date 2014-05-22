function regionUnary = getRegionScoreFromProbImage(...
    neuronProbabilityImage,mitochondriaProbabilityImage,...
    useMitochondriaDetection,addMarginThickness,addMarginValue,...
    setOfCells,sizeR,sizeC,...
    wsIndsForCells,ws,showImg)

% Inputs:
%   neuronProbabilityImage
%   mitochondriaProbabilityImage

disp('Using precomputed probability maps to get region unaries')

neuronProbabilityMap = double(imread(neuronProbabilityImage));
neuronProbabilityMap = neuronProbabilityMap./(max(max(neuronProbabilityMap)));

if(showImg)
    figure;imshow(neuronProbabilityMap);title('Neuron probability map')
end

if(useMitochondriaDetection)
    % overlay mitochondria probabilities (1-P) 
    disp('Mitochondria detection enabled')
    mitochondriaProbabilityMap = double(imread(mitochondriaProbabilityImage));
    mitochondriaProbabilityMap = mitochondriaProbabilityMap./(max(max(mitochondriaProbabilityMap)));
    
    mitoProbInverted = (mitochondriaProbabilityMap -1) * (-1);
    
    % overlay mitochondria as membrane type probability on the neuron
    % probability map
%     neuronProbabilityMap(mitochondriaProbabilityMap>0)...
%         = mitoProbInverted(mitochondriaProbabilityMap>0.2);    

    neuronProbabilityMap = min(neuronProbabilityMap,mitoProbInverted);
    if(showImg)
        figure;imshow(neuronProbabilityMap);
        title('Neuron probability map with mitochondria overlaid')
        
        figure;imshow(mitochondriaProbabilityMap);
        title('Mitochondria probability map')
    end
    
end

% add border to be compatible with filtered image
if(addMarginThickness>0)
    probabilityMapWithMargin = addThickBorder...
                (neuronProbabilityMap,addMarginThickness,addMarginValue);
else
    probabilityMapWithMargin = neuronProbabilityMap;
end


regionUnary = getCellPriors_probability(probabilityMapWithMargin,setOfCells,...
    sizeR,sizeC,wsIndsForCells,ws,showImg);