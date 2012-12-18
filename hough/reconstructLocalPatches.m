function localPatches = reconstructLocalPatches(patchLines,bb)
% Inputs:
% patchLines - a cell array containing the lines (struct) for each patch
% bb - patch size

% Output:
% a cell array containing a matrix corresponding to each image patch whose
% elements are set to 1 according to the detected line segments

% Parameters
stepSize = 0.2;

[rows cols] = size(patchLines);

localPatches = cell(rows,cols);
%x = 1:bb*resolution;

for i=1:rows
    for j=1:cols
        patch = zeros(bb,bb); % removed the resolution setting
        for k = 1:length(patchLines{i,j}) % for each line in this patch
            % create patch
            % TODO: check the existence of point1 with isfield
            if(isfield(patchLines{i,j}(k),'point1'))
                startPt = patchLines{i,j}(k).point1;
                stopPt = patchLines{i,j}(k).point2;
                % x = patchLines{i,j}.point1(1):patchLines{i,j}.point2(1);
                x = startPt(1):stepSize:stopPt(1);
                % x = int8(x);
            else
                break;
            end
            theta = patchLines{i,j}(k).theta;
            rho = patchLines{i,j}(k).rho;
            if(startPt(1)==stopPt(1))
                % if it's a vertical line
                y = startPt(2):stepSize:stopPt(2);
            else
                y = (-cosd(theta)/sind(theta))*x ...
                    + rho/sind(theta);

            end
            x = int8(x);
            y = int8(y);
            y = y + 1;
            biggerInd = find(y>bb);
            y(biggerInd) = bb;

            % draw lines from patchLines{i,j}
            if(size(x,2)>0)                
                patch2 = drawLineInMat(patch,x,y);
                patch = patch + patch2;
            end
        end
        localPatches{i,j} = patch;
    end
end