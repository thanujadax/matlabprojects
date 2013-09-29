% test B-spline

% create set of points (circle)
% center: (5,5), radius 4
a = 5;
b = 5;
r = 4;

x1  = a-r:0.25:a+r;
x2 = fliplr(x1);
y1 = sqrt(r^2 - (x1-a).^2) + b + 1*randn(1);
y2 = -sqrt(r^2 - (x2-a).^2) + b +1*randn(1);

% data vector
data = zeros(numel(x1)*2,2);
data(:,1) = [x1';x2'];
data(1:numel(y1),2) = y1;
data((numel(y1)+1):size(data,1),2) = y2;

% append a copy of the top 4 points at the bottom of the data list
top4 = data(1:4,:);
data = [data; top4];


x = data(:,1);
y = data(:,2);

figure;plot(data(:,1),data(:,2),'*')

% approximate by B-splines
t = linspace( 0, 1, (numel(x)+4) ); % define the parameter t
dt = t(2) - t(1);

% B-spline example 2
% Evaluate and plot a cubic B-spline parametric curve; relies
% on posted Bspline.m routine.

% Evenly spaced knots, with 4 more knots than control points
%t = [0 1 2 3 4 5 6 7 8 9 10 11 12 13];
%vx =    [5 8 9 7 5 2 1 2  5  8];
%vy =    [0 1 3 3 6 4 1 1  3  0];

% evaluate at a number of points for plotting
qx = Bspline(t, x, dt);
qy = Bspline(t, y, dt);

% plot(x, y, '--', qx, qy, '-');
% xlabel('t');
% ylabel('y');
% legend('control polygon', 'B-spline curve');
figure;plot(qx,qy)
%% spline curve by uniform subdivision
order = 3;
values = spcrv(data',order); 
figure;
plot(data(:,1),data(:,2),'*')
hold on, plot(values(1,:),values(2,:)), hold off
