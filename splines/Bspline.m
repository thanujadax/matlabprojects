function q = Bspline(t, v, dt)
% Evaluate the piecewise cubic B-spline curve at intervals of dt
% and return the resulting vector of points
%  has knots t(i) and control points v(i).
%  IMPORTANT: 
%   (1) there should be 4 more knots t(i) than points v(i).
%   (2) This function evaluates the B-spline curve over [t_4, t_{m-3}],
%       where m is the number of knots. (Note the spline is not 
%       well-defined outside this interval.)


m = size(t, 2);  % number of knots
i = 4;           % index of first knot interval over which curve is evaluated
ptNum = 0;       % counter for points generated

% Loop through all evaluation values
for u = t(4) : dt: t(m-3);
   ptNum = ptNum + 1;
   % check if u value has moved to the next knot interval
   % include small tolerance on knots to avoid round-off error in comparisons.
   while (u > t(i+1)+1.0e-10)
       i = i + 1;
   end

   % Now evaluate the spline at u using the deBoor algorithm.
   % Start with the relevant control points.
   % w used here to simplify indices.
   w = i - 4;
   for j = 1 : 4
     qq(j) = v(w + j);
   end
   for j = 1 : 3
     for k = 1 : 4 - j
        qq(k) = (t(w + k + 4) - u)/(t(w + k + 4) - t(w + k + j)) * qq(k) + ...
             (u - t(w + k + j))/(t(w + k + 4) - t(w + k + j)) * qq(k+1);
     end
   end
   % Create vector of points on the B-spline curve.
   q(ptNum) = qq(1);
end