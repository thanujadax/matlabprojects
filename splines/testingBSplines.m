% test B-spline

% create set of points (circle)
% center: (5,5), radius 4
a = 5;
b = 5;
r = 4;

x  = a-r:a+r;
y1 = sqrt(r^2 - (x-a).^2) + b;
y2 = -sqrt(r^2 - (x-a).^2) + b;

figure;
plot(x,y1)
hold on
plot(x,y2)
title('points joined with straight lines')
% approximate by B-splines

