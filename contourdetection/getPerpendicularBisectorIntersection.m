function [x,y] = getPerpendicularBisectorIntersection(edge1_pixels,edge2_pixels,sizeR,sizeC,...
                    edges2nodes,nodeInds)

% get the midpoint
midpoint1 = getMidpoint(edge1_pixels);
midpoint2 = getMidpoint(edge2_pixels);


% get the tangent at the midpoint



% draw perpendiculars to the tangents at the midpoints
% get the intersection

end

function midpoint = getMidpoint(edgepixels)

numpixels = numel(edgepixels);
midpoint_pos = ceil((numpixels+1)/2);
midpoint = edgepixels(midpoint_pos);

end

function [mdegrees,c] = getTangent(edgepixels,midpoint,edges2nodes,edgeInds)
    numEdgePix = numel(edgepixels)
    if(numEdgePix>2)
        midpos = find(edgepixels==midpoint);
        pt1 = edgepixels(midpos-1);
        pt2 = edgepixels(midpos+1);    
    else
        % there're only 2 or 1 pixels for this edge. then we need the nodepixels
        % on either side
        pt1 = edges2nodes()
    end
    [x1,y1] = ind2sub([sizeR,sizeC],pt1);
    [x2,y2] = ind2sub([sizeR,sizeC],pt2);
    [mdegrees,c] = getStraightLineEqn(x1,y1,x2,y2); 
    
end

function [mdegrees,c] = getStraightLineEqn(x1,y1,x2,y2)
    mdegrees = atan2d((y2-y1),(x2-x1));
    if(mdegrees==90)
        c = 0; % x has just one value for all values of y
    else
        c = y1 - x1 .* tand(mdegrees);
    end
end