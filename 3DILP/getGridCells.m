function [griIDs_sectionIDs_rootPixIDsRel,numEdgesY,numEdgesX]...
            = getGridCells...
            (numR,numC,numZ,gridResX,gridResY,gridResZ)

% Input:
%   inputImageStackMat - 3D matrix of the image stack. type: uint8
%       includes membrane cells around the borders

% Output:
%   griIDs_sectionIDs_rootPixIDs - matrix where each col is suggested by name
%   rootPixID is the pixInd of the start pixel (0,0,0) of each cell,
%   wrt each slice (relative coordinates in each slice)


% grid dimensions (numCells in each dim)

numEdgesX = floor((numC-1)/(gridResX));
marginPix_X = mod((numC-1),(gridResX));
gridStartX = floor(marginPix_X/2);
gridStartX = max(1,gridStartX);

numEdgesY = floor((numR-1)/(gridResY));
marginPix_Y = mod((numC-1),(gridResY));
gridStartY = floor(marginPix_Y/2);
gridStartY = max(1,gridStartY);

% define node positions
x_pos = gridStartX:gridResX:numC;
y_pos = gridStartY:gridResY:numR;

% Here, a grid node corresponds to the root (0,0,0) of each cube.
% Therefore, we get rid of the final 'grid node' which is the one at the
% boundary
x_pos(end) = [];
y_pos(end) = [];

[nodeX, nodeY] = meshgrid(x_pos,y_pos);

nodePixList_section = sub2ind([numR numC],nodeY,nodeX);
numGridCellsPerSection = numEdgesX * numEdgesY;
totGridCells = numGridCellsPerSection * numZ;

griIDs_sectionIDs_rootPixIDsRel = uint8(zeros(totGridCells,3));
gridCounter = 0;
for k=1:numZ
    for j=1:numEdgesX
        for i=1:numEdgesY
            gridCounter = gridCounter + 1;
            % gridCellInd
            griIDs_sectionIDs_rootPixIDsRel(gridCounter,1) = gridCounter;
            % section_ID
            griIDs_sectionIDs_rootPixIDsRel(gridCounter,2) = k;
            % rootPixInd
            griIDs_sectionIDs_rootPixIDsRel(gridCounter,3)...
                = sub2ind([numR numC],nodeY(i),nodeX(j));
        end
    end
end