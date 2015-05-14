classdef imageObj < hgsetget
    
    properties (SetAccess = protected)
        imageMat
        
    end
    
    methods
        % constructor
        function obj = imageObj(numRows,numCols)
            obj.imageMat = zeros(numRows,numCols);
        end
        
        function changeImagePixel(obj,coordinates)
            obj.imageMat(ceil(coordinates(2)),ceil(coordinates(1))) = 1;                        
        end
       
    end
    
end