function [x,y] = getPerpendicularBisectorIntersection(edge1_pixels,edge2_pixels,...
    edge1_id,edge2_id, sizeR,sizeC,edges2nodes,nodeInds,edgeIndsAll)

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
normal1 = getNormalToLine(tangent1,midxy1);
normal2 = getNormalToLine(tangent2,midxy2);
% get the intersection
[x,y] = getIntersectionOfTwoLines(normal1,normal2);
x = floor(x);
y = floor(y);
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
        edgeNodes_listInds = edges2nodes(edgeIndsAll==edgeID,:);
        
        pt1 = nodeInds(edgeNodes_listInds(1));
        pt2 = nodeInds(edgeNodes_listInds(2));        
    end
    pt3 = edgepixels(midpos);
    
    % create structure array to contain the 3 points
    p = struct([]);
    [p(1).y, p(1).x] = ind2sub([sizeR,sizeC],pt1);
    [p(2).y, p(2).x] = ind2sub([sizeR,sizeC],pt2);
    [p(3).y, p(3).x] = ind2sub([sizeR,sizeC],pt3);
    
    lineEqn = getStraightLineEqn(p); 
    
end

function lineEqn = getStraightLineEqn(p)
    % p contains 3 points. but we need just 2.
    syms a b c x y eq1 eq2 eq3 eq
    eq = a*x + b*y + c;
    % let b = 1;
    eq1 = subs(eq,[x,y,b],[p(1).x,p(1).y,1]);
    eq2 = subs(eq,[x,y,b],[p(2).x,p(2).y,1]);
    Eqn = solve(eq1,eq2,a,c);
    if(isempty(Eqn))
        % then let b = 0, a=1
        eq3 = subs(eq,[x,y,b,a],[p(1).x,p(1).y,0,1]);
        Eqn = solve(eq3,c);
        lineEqn.c = Eqn;
        lineEqn.a = 1;
        lineEqn.b = 0;
    else
        lineEqn.a = Eqn.a;
        lineEqn.b = 1;
        lineEqn.c = Eqn.c;
    end
end

function normalLine = getNormalToLine(lineEqn,midxy)
    syms a b c a0 b0 c0 x y eq1 eq2
    eq1 = a*a0 + b*b0;
    eq2 = a*x + b*y + c;
    eq3 = subs(eq1,[a0,b0,b],[lineEqn.a,lineEqn.b,1]);
    eq4 = subs(eq2,[x,y,b],[midxy.x,midxy.y,1]); % subs b = 1.
    % then a = -gradient, c = -intercept
    nline = solve(eq3,eq4,a,c);
    if(isempty(nline))
        eq4 = subs(eq2,[x,y,a,b],[midxy.x,midxy.y,1,0]); % subs a=1,b=0
        % then the eqn is in the form x = c
        normalLine.a = 1;
        normalLine.b = 0;
        nline = solve(eq4,c);
        normalLine.c = nline;
    else
        normalLine.a = nline.a;
        normalLine.b = 1;
        normalLine.c = nline.c;
    end
end

function [x,y] = getIntersectionOfTwoLines(line1,line2)
% normalize lines
syms a b c x y eq1 eq2 eq
eq = a*x + b*y + c;
eq1 = subs(eq,[a,b,c],[line1.a,line1.b,line1.c]);
eq2 = subs(eq,[a,b,c],[line2.a,line2.b,line2.c]);
pt = solve(eq1,eq2,x,y);
% TODO: convert sym to double
if(~isempty(pt))
    x = double(pt.x);
    y = double(pt.y);
else
   x = [];
   y = [];
end
end