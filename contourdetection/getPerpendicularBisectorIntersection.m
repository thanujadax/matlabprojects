function [x,y] = getPerpendicularBisectorIntersection(edge1_pixels,edge2_pixels,...
    edge1_id,edge2_id, sizeR,sizeC,edges2nodes,nodeInds)

x = [];
y = [];

% get the midpoint
midpoint1 = getMidpoint(edge1_pixels);
midpoint2 = getMidpoint(edge2_pixels);
% convert to xy coordinates
[midxy1.y,midxy1.x] = ind2sub([sizeR sizeC],midpoint1);
[midxy2.y,midxy2.x] = ind2sub([sizeR sizeC],midpoint2);

% get the tangent at the midpoint
% (eqn of the lines)
tangent1 = getTangentAtPoint(edge1_pixels,midpoint1,edges2nodes,edgeIndsAll,...
    sizeR,sizeC,edge1_id,nodeInds);
tangent2 = getTangentAtPoint(edge2_pixels,midpoint2,edges2nodes,edgeIndsAll,...
    sizeR,sizeC,edge2_id,nodeInds);

% draw perpendiculars to the tangents at the midpoints
normal1 = getNormalToLine(line1,midxy1);
normal2 = getNormalToLine(line2,midxy2);
% get the intersection
[x,y] = getIntersectionOfTwoLines(normal1,normal2);
end

function midpoint = getMidpoint(edgepixels)

numpixels = numel(edgepixels);
midpoint_pos = ceil((numpixels+1)/2);
midpoint = edgepixels(midpoint_pos);

end

function lineEqn = getTangentAtPoint(edgepixels,midpoint,edges2nodes,edgeIndsAll,...
    sizeR,sizeC,edgeID,nodeInds)

    midpos = find(edgepixels==midpoint);
    numEdgePix = numel(edgepixels);
    if(numEdgePix>2)
        pt1 = edgepixels(midpos-1);
        pt2 = edgepixels(midpos+1);        
    else
        % there're only 2 or 1 pixels for this edge. then we need the nodepixels
        % on either side
        edgeNodes_listInds = edges2nodes(edgeIndsAll==edgeID);
        pt1 = nodeInds(edgeNodes_listInds(1));
        pt2 = nodeInds(edgeNodes_listInds(2));        
    end
    pt3 = edgespixels(midpos);
    
    % create structure array to contain the 3 points
    p = struct([]);
    [p(1).y, p(1).x] = ind2sub([sizeR,sizeC],pt1);
    [p(2).y, p(2).x] = ind2sub([sizeR,sizeC],pt2);
    [p(3).y, p(3).x] = ind2sub([sizeR,sizeC],pt3);
    
    lineEqn = getStraightLineEqn(p); 
    
end

function lineEqn = getStraightLineEqn(p)
    % p contains 3 points. solve ax + by + c = 0 to find a, b, c.
    syms a b c x y eq1 eq2 eq3 eq
    eq = a*x + b*y + c;
    eq1 = subs(eq,[x,y],[p(1).x,p(1).y]);
    eq2 = subs(eq,[x,y],[p(2).x,p(2).y]);
    eq3 = subs(eq,[x,y],[p(3).x,p(3).y]);
    [lineEqn.a,lineEqn.b,lineEqn.c] = solve(eq1,eq2,eq3,a,b,c);
end

function normalLine = getNormalToLine(lineEqn,midxy)
    syms a b c a0 b0 c0 x y eq1 eq2
    eq1 = a*a0 + b*b0;
    eq2 = a*x + b*y + c;
    eq3 = subs(eq1,[a0,b0,b],[lineEqn.a,lineEqn.b,1]);
    eq4 = subs(eq2,[x,y,b],[midxy.x,midxy.y,1]); % subs b = 1.
    % then a = -gradient, c = -intercept
    [a,c] = solve(eq3,eq4,a,c);
    if(a==inf || a==(-inf))
        eq4 = subs(eq2,[x,y,a,b],[midxy.x,midxy.y,1,0]); % subs a=1,b=0
        % then the eqn is in the form x = c
        normalLine.a = 1;
        normalLine.b = 0;
        normalLine.c = solve(eq4,c);
    else
        normalLine.a = a;
        normalLine.b = 1;
        normalLine.c = c;
    end
end

function [x,y] = getIntersectionOfTwoLines(line1,line2)
% normalize lines
syms a b c x y eq1 eq2 eq
eq = a*x + b*y + c;
eq1 = subs(eq,[a,b,c],[line1.a,line1.b,line1.c]);
eq2 = subs(eq,[a,b,c],[line2.a,line2.b,line2.c]);
[x,y] = solve(eq1,eq2,x,y);
end