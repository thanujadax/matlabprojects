function lines = patchLinesToGlobal(patchLines,slidingDist,bb)

% TODO: correct the shifting for the sliding distance!!!

% Inputs:
% patchLines - 2D cell array. each cell contains an array of structs each
% containing the information about local line segments
%   point1 = [x1 y1]
%   point2 = [x2 y2]
%   theta = angle in degrees
%   rho = position

% translate each patchLine according to the patch position, patch size and the sliding
% distance.

% output:
% imOut - reconstruction


[rows cols] = size(patchLines);
m = 1;

for i=1:rows
    for j=1:cols
        for k = 1:length(patchLines{i,j}) % for each line in this patch
            if(isfield(patchLines{i,j},'point1'))
                % non-empty struct found
                % shifting (x,y)
                xshift = (j-1)*(bb-slidingDist);
                yshift = (i-1)*(bb-slidingDist);
                
                point1 = patchLines{i,j}(k).point1;
                point1 = point1 + [xshift yshift];
                lines(m).point1 = point1;
                
  
                point2 = patchLines{i,j}(k).point2;
                point2 = point2 + [xshift yshift];
                lines(m).point2 = point2;
                               
                % shift rho
                rho = patchLines{i,j}(k).rho;
                theta = patchLines{i,j}(k).theta;
                xr = rho*cosd(theta) + (j-1)*(bb-slidingDist);
                yr = rho*sind(theta) + (i-1)*(bb-slidingDist);
                newRho = sqrt(xr^2 + yr^2);

                lines(m).rho = newRho;
                lines(m).theta = theta;
                m = m + 1;          
            else
                % empty struct. go to next iteration
                break;
            end
            
        end
    end
end


